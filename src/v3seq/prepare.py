#!/usr/bin/env python3

import sys
import os
import shlex
import logging
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pkg_resources import resource_filename

cons_file = resource_filename(__name__, 'db/cons.faa')

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    from itertools import zip_longest
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def filter_reads(filename, max_n, min_len=129):
    """Use seqtk and Biopython to trim and filter low quality reads."""

    # run seqtk trimfq to trim low quality ends
    logging.info('Trimming reads with seqtk')
    r1 = 'seqtk trimfq %s | seqtk sample - %d | seqtk seq -L %d -A - > high_quality.fasta' % (filename, max_n, min_len)
    subprocess.call(r1, shell=True, universal_newlines=True)
    return 'high_quality.fasta'


def blast_reads(read_group):
    """Use blastx to align reads to V3 consensus"""
    SeqIO.write(read_group, 'tmp.fasta', 'fasta')
    cml = "blastx -task blastx-fast -query tmp.fasta -subject %s \
    -out out.tsv -num_alignments 1 -evalue 1E-4 \
    -outfmt '6 qseqid sseqid pident qcovs score length mismatch gapopen qstart qend sstart send'" % (cons_file)
    subprocess.call(cml, shell=True)
    os.remove('tmp.fasta')
    als = pd.read_table('out.tsv', header=None,
                        names=['qseqid', 'sseqid', 'pident', 'qcovs', 'score', 'length', 'mismatch', 'gapopen',
                               'qstart', 'qend', 'sstart', 'send'])
    os.remove('out.tsv')
    covering = als[als.send - als.sstart > 33]
    if covering.empty:
        return []
    covering['qseqid'] = covering.apply(
        lambda x: x['qseqid'] if x['qstart'] < x['qend'] else str(x['qseqid']) + ':REV', axis=1)
    return covering.qseqid.tolist()

def main(filein):
    """What the main does."""
    logging.info('obtain high quality reads')
    hq = filter_reads(filein, 10000)
    hq_reads = SeqIO.parse(hq, 'fasta')
    covering_reads = set([])
    logging.info('blast reads in batches until enough are found')
    for i, group in enumerate(grouper(1000, hq_reads)):
        print('blast call %d' % i)
        _ = blast_reads(group)
        covering_reads.update(_)
        if len(covering_reads) >= 50:
            break

    cr = [s for s in SeqIO.parse(hq, 'fasta') if s.id in covering_reads]
    SeqIO.write(cr, 'out.fasta', 'fasta')
    crev = [SeqRecord(s.seq.reverse_complement(), s.id, description='reversed') for s in SeqIO.parse(hq, 'fasta') \
            if '%s:REV' % s.id in covering_reads]
    SeqIO.write(crev, 'revout.fasta', 'fasta')
    subprocess.call('cat out.fasta revout.fasta > v3reads.fasta', shell=True)
    cml = shlex.split('muscle -in v3reads.fasta -out msa.fasta -quiet')
    subprocess.call(cml)
    cml = shlex.split('cons -sequence msa.fasta -outseq v3cons.fasta -plurality 0.9 -identity 10 -name v3cons -auto')
    subprocess.call(cml)
    for f in ['high_quality.fasta', 'out.fasta', 'revout.fasta']:
        os.remove(f)

if __name__ == '__main__':
    main(sys.argv[1])
