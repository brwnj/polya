#!/usr/bin/env python
# encoding: utf-8
"""
Since counts are in the form of bedgraph, counts as a result of this script
are unstranded. Any input should be strand specific.
"""
import subprocess as sp
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def map_counts(bedgraph, bed):
    """print counts to stdout using bedtools map"""
    # overlap and get counts
    # add appropriate gene name
    # add indeces for sorting
    # sort by gene name then site number
    # output dexseq columns
    cmd = """bedtools map -c 4 -o max -null 0 -a %s -b %s |\
              awk '{split($4, symbol, "|"); split(symbol[1], gene, ".");print gene[3]":"$4"\t"$7}' |\
              awk '{split($1, full, ":"); split(full[2], site, "."); print full[1]"\t"site[4]"\t"$0}' |\
              sort -k1,1 -k2,2n |\
              awk '{print $3"\t"$4}'"""\
              % (bed, bedgraph)
    
    sp.call(cmd, shell=True)

def main(counts, sites):
    map_counts(counts, sites)

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("counts", help="read counts for a sample as bedgraph")
    p.add_argument("sites", help="merged polya sites as bed")
    args = p.parse_args()
    main(args.counts, args.sites)
