#!/usr/bin/env python3
"""Everything happens here."""


import sys
import os
import shlex
import logging
import subprocess
import re
from collections import Counter

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
    """Take a MSA in fasta format where headers are like read_10-count_15 and return:
       - a data frame with position, nucleotide, frequencies;
       - a Counter with haplotypes."""
    msa = AlignIO.read(filename, 'fasta')
    n = msa.get_alignment_length()  # columns, number of positions
    pos = []
    nt = []
    freq = []
    for j in range(n):  # iterate over columns
        tmp_nt = []
        for i, seq in enumerate(msa):
            read_count = int(seq.id.split('_')[-1])  # extract count from seq.id
            tmp_nt.extend([msa[i, j]] * read_count)  # add as many nt as needed
        for b, counts in Counter(tmp_nt).items():
            pos.append(j)
            nt.append(b)
            freq.append(float(counts) / len(tmp_nt))
    df = pd.DataFrame({'pos': pos, 'nt': nt, 'freq': freq})

    all_reads = []
    for s in msa:
        read_count = int(s.id.split('_')[-1])
        all_reads.extend([str(s.seq)] * read_count)

    haplos = Counter(all_reads)

    return df, haplos, len(tmp_nt)


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


def remove_matching_reads(filename, cont_file):
    """Remove reads with matches anything in cont_file."""
    if not os.path.exists(cont_file + '.bwt'):
        cml = shlex.split('bwa index %s' % cont_file)
        subprocess.call(cml)
    cml = 'bwa mem -t 2 %s %s 2> /dev/null | samtools view -f 4 -h - | samtools bam2fq - ' % (cont_file, filename)
    cml += '| seqtk seq -A - > clean_reads.fasta'

    subprocess.call(cml, shell=True)
    return 'clean_reads.fasta'


def filter_reads(filename, max_n=100000, min_len=129):
    """Use seqtk to trim, subsample, filter short reads."""
    # run seqtk trimfq to trim low quality ends
    logging.info('Trimming reads with seqtk, subsample, and delete reads shorter than %d', min_len)
    r1 = 'seqtk trimfq %s | seqtk seq -L %d | seqtk sample - %d > high_quality.fastq' % (filename, min_len, max_n)
    subprocess.call(r1, shell=True, universal_newlines=True)
    return 'high_quality.fastq'


def blast_reads(read_group):
    """Use blastx to align reads to V3 consensus."""
    try:
        tot_seqs = SeqIO.write(read_group, 'tmp.fasta', 'fasta')
    except AttributeError:  # captures when last read_group is filled with None at the end
        return []
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
    covering = als.copy()
    covering = covering[(covering.send - covering.sstart > 33)]
    if covering.empty:
        return []
    value = covering.apply(
        lambda x: '%s:FWD:%d:%d' % (x['qseqid'], x['qstart'], x['qend']) if x['qstart'] < x['qend'] else
        '%s:REV:%d:%d' % (x['qseqid'], x['qend'], x['qstart']), axis=1)
    covering.loc[:, 'qseqid'] = value
    logging.info('Found %d covering reads', covering.shape[0])
    return covering.qseqid.tolist()


def extract_reads(reads_list, reads_file):
    """Extract reads from reads_file according to id:orientation:start:end in reads_list."""
    all_reads = SeqIO.to_dict(SeqIO.parse(reads_file, 'fasta'))
    save_reads = []
    n1 = n2 = 0
    for read_info in reads_list:
        sid, orientation, start, end = re.search(r'(.*):([A-Z]{3}):(\d*):(\d*)$', read_info).group(1, 2, 3, 4)
        if orientation == 'FWD':
            save_reads.append(str(all_reads[sid][int(start) - 1:int(end)].seq))
            n1 += 1
        elif orientation == 'REV':
            rh = str(all_reads[sid][int(start) - 1:int(end)].seq.reverse_complement())
            save_reads.append(rh)
            n2 += 1
    count_reads = Counter(save_reads)
    sr = [SeqRecord(Seq(read), id='read_%d-count_%d' % (i, count), description='')
          for i, (read, count) in enumerate(count_reads.items())]
    return sr, n1, n2


def main(filein, min_reads=150, n_group=2000):
    """What the main does."""
    from random import sample
    assert os.path.exists(filein)
    hq = filter_reads(filein)
    logging.info('remove matching reads')
    no_pol = remove_matching_reads(hq, cont_file)
    # no_pol = 'clean_reads.fasta'
    no_pol_reads = list(SeqIO.parse(no_pol, 'fasta'))
    no_pol_reads = sample(no_pol_reads, k=len(no_pol_reads))
    covering_reads = set([])
    logging.info('blast reads in batches until enough are found')
    total_blasted = 0
    for i, group in enumerate(grouper(n_group, no_pol_reads)):
        if i > 2 and len(covering_reads) < 20:
            sys.exit('not enough reads covering V3 were found')
        logging.info('blast call %d', i + 1)
        _ = blast_reads(group)
        covering_reads.update(_)
        total_blasted += n_group
        logging.info('this blast: %d covering  out of %d total - %3.2f %%', len(_), n_group,
                     100 * float(len(_)) / n_group)
        logging.info('cumulative: %d covering  out of %d total - %3.2f %%', len(covering_reads), total_blasted,
                     100 * float(len(covering_reads)) / total_blasted)
        if len(covering_reads) >= min_reads:
            break

    logging.info('covering_reads used in MSA: %d out of %d blasted (%3.2f %%)', len(covering_reads), total_blasted,
                 100 * float(len(covering_reads)) / total_blasted)
    cov_reads, n_fwd, n_rev = extract_reads(covering_reads, no_pol)

    SeqIO.write(cov_reads, 'v3reads.fasta', 'fasta')
    logging.info('%d covering reads in forward orientation', n_fwd)
    logging.info('%d covering reads in reverse orientation', n_rev)
    if n_fwd + n_rev < min_reads:
        logging.error('Not enough reads: %d', n_fwd + n_rev)
        sys.exit('Not enough reads: %d' % (n_fwd + n_rev))

    no_singleton_reads = [s for s in SeqIO.parse('v3reads.fasta', 'fasta') if int(s.id.split('_')[-1]) > 1]
    SeqIO.write(no_singleton_reads, 'v3reads_no_singleton.fasta', 'fasta')

    cml = shlex.split('muscle -in v3reads_no_singleton.fasta -out msa.fasta -quiet')
    subprocess.call(cml)

    df, haplotypes, support = msa_2_df('msa.fasta')
    logging.info('Haplotypes supported by %d reads out of %d: %3.1f%%',
                 support, n_fwd + n_rev, 100.0 * support / (n_fwd + n_rev))
    cons_seq = df_2_ambiguous_sequence(df)
    SeqIO.write([SeqRecord(Seq(cons_seq), id='v3_consensus', description='')], 'v3cons.fasta', 'fasta')

    haps = []
    hi = 1  # counter for haplotypes, used in fasta file
    accounted_f = 0.0  # keep track of the cumulative accounted frequency
    tot_reads = sum(haplotypes.values())
    for h, support in haplotypes.most_common():
        f = round(float(support) / tot_reads, 2)
        accounted_f += f
        sr = SeqRecord(Seq(h), id='v3_haplotype_%d-support_%3.2f' % (hi, f), description='')
        haps.append(sr)
        hi += 1

    SeqIO.write(haps, 'v3haplotypes.fasta', 'fasta')
    for f in ['high_quality.fastq', 'clean_reads.fasta']:
        os.remove(f)
    logging.info('Haplotypes written to haplotypes.fasta')


if __name__ == '__main__':
    main(sys.argv[1])
