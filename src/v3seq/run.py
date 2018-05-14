#!/usr/bin/env python3
"""Everything happens here."""


import sys
import os
import shlex
import logging
import subprocess
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename

d2a = {'AG': 'R', 'CT': 'Y', 'AC': 'M', 'GT': 'K', 'CG': 'S', 'AT': 'W',
       'ACT': 'H', 'CGT': 'B', 'ACG': 'V', 'AGT': 'D', 'ACGT': 'N'}

cons_file = resource_filename(__name__, 'db/v3cons.faa')
cont_file = resource_filename(__name__, 'db/pol_sequences.fasta')

try:
    n_proc = min(os.cpu_count(), 8)
except NotImplementedError:
    n_proc = 2
logging.info('%d cores that will be used', n_proc)

def grouper(n, iterable, fillvalue=None):
    """Group an iteranle n elements at a time, from itertools examples.

    grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    """
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def msa_2_df(filename):
    """Take a MSA in fasta format and return:
       - a data frame with position, nucleotide, frequencies;
       - a Counter with haplotypes."""
    from collections import Counter
    msa = AlignIO.read(filename, 'fasta')
    m = len(msa)  # rows, number of sequences
    n = msa.get_alignment_length()  # columns, number of positions
    j1 = 0
    j2 = n
    for j in range(n):
        c = Counter(msa[:, j])
        if c['-'] < 5:
            j1 = j
            break
    for j in range(n - 1, -1, -1):
        c = Counter(msa[:, j])
        if c['-'] < 5:
            j2 = j
            break
    window = msa[:, j1:j2 + 1]
    haplos = Counter((str(window[i, :].seq) for i in range(m)))

    pos = []
    nt = []
    freq = []
    for j in range(j1, j2):
        for b, counts in Counter(msa[:, j]).items():
            pos.append(j)
            nt.append(b)
            freq.append(float(counts) / m)
    df = pd.DataFrame({'pos': pos, 'nt': nt, 'freq': freq})
    return df, haplos


def df_2_ambiguous_sequence(df_in):  # , cov_df=None):
    """Take a DataFrame with positions, nucleotides and frequencies and returns a sequence.

    If the frequency is above a certain threshold, write the sequence with wobble bases.
    """
    assert 'freq' in df_in.columns
    # select calls with freq > 15%
    df_in = df_in[df_in['freq'] >= 0.15]
    # aggregate calls for the same position
    all_nt = df_in.groupby(['pos']).agg({'nt': lambda x: ''.join(sorted(x))})
    # create a columng of ambiguous bases
    value = all_nt.apply(lambda row: d2a.get(row['nt'], row['nt'][0]), axis=1)
    all_nt.loc[:, 'ambi'] = value
    all_nt.reset_index(inplace=True)
    # if not cov_df is None:
    #     full_df = pd.merge(all_nt, cov_df, on='pos', how='left')
    #     full_df.loc[full_df['coverage'] < coverage_threshold, 'ambi'] = 'N'
    return ''.join(all_nt.ambi.tolist())


def remove_matching_reads(filename, contaminant_file):
    """Remove reads with matches anything in contaminant_file."""
    if not os.path.exists(cont_file + '.bwt'):
        cml = shlex.split('bwa index %s' % cont_file)
        subprocess.call(cml)
    cml = 'bwa mem -t 2 %s %s | samtools view -f 4 -h - | samtools bam2fq - | seqtk seq -A - > clean_reads.fasta' % \
        (contaminant_file, filename)
    subprocess.call(cml, shell=True)
    return 'clean_reads.fasta'


def filter_reads(filename, max_n=100000, min_len=129):
    """Use seqtk and Biopython to trim and filter low quality reads."""
    # run seqtk trimfq to trim low quality ends
    logging.info('Trimming reads with seqtk')
    r1 = 'seqtk trimfq %s | seqtk sample - %d | seqtk seq -L %d - > high_quality.fastq' % (filename, max_n, min_len)
    subprocess.call(r1, shell=True, universal_newlines=True)
    return 'high_quality.fastq'


def blast_reads(read_group):
    """Use blastx to align reads to V3 consensus."""
    tot_seqs = SeqIO.write(read_group, 'tmp.fasta', 'fasta')

    max_n = (tot_seqs / n_proc) + 1
    # We want to split in n_proc processors, so each file has at most
    # (tot_seqs / n_proc) + 1 reads
    cml = "awk -v \"MAX_N=%d\" \'BEGIN {n_seq=0;} /^>/ \
    {if(n_seq %% %d == 0){file=sprintf(\"splitted_clean_%%d.fasta\", n_seq/%d);} \
    print >> file; n_seq++; next;} { print >> file; }' %s" % (max_n, max_n, max_n, 'tmp.fasta')
    subprocess.call(cml, shell=True)

    if sys.platform.startswith('linux'):
        xargs_thread = 0  # means on all available cores, caution
    elif sys.platform.startswith('darwin'):
        xargs_thread = n_proc  # darwin xargs does not accept -P 0
    else:
        logging.debug('could not detect system platform: runnning on %d cores', n_proc)
        xargs_thread = n_proc

    cml = 'seq 0 %s | xargs -P %d -I {} blastx -task blastx-fast -subject %s \
           -query splitted_clean_{}.fasta -out tmp_{}.tsv -num_alignments 1 -evalue 1E-4 \
           -outfmt \'6 qseqid sseqid pident qcovs score length mismatch gapopen qstart qend sstart send\'' \
        % (n_proc - 1, xargs_thread, cons_file)
    logging.debug('running blast in parallel now')
    subprocess.call(cml, shell=True)

    subprocess.call('rm tmp.fasta splitted_clean_*fasta', shell=True)
    subprocess.call('cat tmp*.tsv > out.tsv', shell=True)
    subprocess.call('rm tmp*.tsv', shell=True)
    als = pd.read_table('out.tsv', header=None,
                        names=['qseqid', 'sseqid', 'pident', 'qcovs', 'score', 'length', 'mismatch', 'gapopen',
                               'qstart', 'qend', 'sstart', 'send'])
    os.remove('out.tsv')
    covering = als[als.send - als.sstart > 33]
    if covering.empty:
        return []
    value = covering.apply(
        lambda x: x['qseqid'] if x['qstart'] < x['qend'] else str(x['qseqid']) + ':REV', axis=1)
    covering.loc[:, 'qseqid'] = value
    return covering.qseqid.tolist()

def main(filein, min_reads=500, n_group=2000):
    """What the main does."""
    assert os.path.exists(filein)
    logging.info('obtain high quality reads')
    hq = filter_reads(filein)
    logging.info('remove matching reads')
    no_pol = remove_matching_reads(hq, cont_file)
    no_pol = 'clean_reads.fasta'
    no_pol_reads = SeqIO.parse(no_pol, 'fasta')
    covering_reads = set([])
    logging.info('blast reads in batches until enough are found')
    total_blasted = 0
    for i, group in enumerate(grouper(n_group, no_pol_reads)):
        logging.info('blast call %d', i + 1)
        _ = blast_reads(group)
        covering_reads.update(_)
        total_blasted += n_group
        logging.debug('%d covering reads out of %d blasted - %3.2f %%', len(_), n_group, 100 * float(len(_)) / n_group)
        logging.debug('covering_reads found so far: %d out of %d total - %3.2f %%', len(covering_reads), total_blasted,
                      100 * float(len(covering_reads)) / total_blasted)
        if len(covering_reads) >= min_reads:
            break

    logging.info('covering_reads used in MSA: %d out of %d blasted (%3.2f %%)', len(covering_reads), total_blasted,
                 100 * float(len(covering_reads)) / total_blasted)
    cr = [s for s in SeqIO.parse(no_pol, 'fasta') if s.id in covering_reads]
    n1 = SeqIO.write(cr, 'out.fasta', 'fasta')
    logging.info('%d covering reads in forward orientation', n1)
    crev = [SeqRecord(s.seq.reverse_complement(), s.id, description='reversed') for s in SeqIO.parse(no_pol, 'fasta') \
            if '%s:REV' % s.id in covering_reads]
    n2 = SeqIO.write(crev, 'revout.fasta', 'fasta')
    logging.info('%d covering reads in reverse orientation', n2)
    if n1 + n2 < min_reads:
        logging.error('Not enough reads: %d', n1 + n2)
        sys.exit('Not enough reads: %d' % (n1 + n2))
    subprocess.call('cat out.fasta revout.fasta > v3reads.fasta', shell=True)
    cml = shlex.split('muscle -in v3reads.fasta -out msa.fasta -quiet')
    subprocess.call(cml)

    df, haplotypes = msa_2_df('msa.fasta')
    cons_seq = df_2_ambiguous_sequence(df)
    SeqIO.write([SeqRecord(Seq(cons_seq), id='v3_consensus', description='')], 'v3cons.fasta', 'fasta')
    haps = []
    hi = 1  # counter for haplotypes, used in fasta file
    accounted_f = 0.0  # keep track of the cumulative accounted frequency
    tot_reads = sum(haplotypes.values())
    for h, support in haplotypes.most_common():
        f = round(float(support) / tot_reads, 2)
        if f < 0.05:
            break
        accounted_f += f
        sr = SeqRecord(Seq(h), id='v3_haplotype_%d-support_%3.2f' % (hi, f), description='')
        haps.append(sr)
        hi += 1
    logging.info('Total frequency of haplotypes below 5%%: %f', 1 - accounted_f)
    SeqIO.write(haps, 'haplotypes.fasta', 'fasta')
    for f in ['high_quality.fastq', 'out.fasta', 'revout.fasta', 'clean_reads.fasta']:
        os.remove(f)
    logging.info('Haplotypes written to haplotypes.fasta')

if __name__ == '__main__':
    main(sys.argv[1])
