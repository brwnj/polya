#!/usr/bin/env python
# encoding: utf-8
"""
Since counts are in the form of bedgraph, counts as a result of this script
are unstranded. Any input should be strand specific.
"""
import os
import tempfile
import subprocess as sp
# from cruzdb import Genome

# def get_sizes(species):
#     tmp = tempfile.mkstemp(suffix=".bed")[1]
#     g = Genome(db=species)
#     df = g.dataframe("chromInfo")
#     df.to_csv(tmp, cols=['chrom','sizes'],
#                 sep="\t", header=False, index=False)
#     return tmp

def add_slop(bed, sizes, n):
    """add slop onto polya sites."""
    tmp = tempfile.mkstemp(suffix=".bed")[1]
    cmd = "bedtools slop -b %d -i %s -g %s |\
            bedtools sort -i - > %s" % (n, bed, sizes, tmp)
    sp.call(cmd, shell=True)
    return tmp

def map_counts(bed, bedgraph):
    """print counts to stdout using bedtools map"""
    cmd = """bedtools map -c 4 -o max -null 0 -a %s -b %s |\
              awk '{split($4, symbol, "|"); split(symbol[1], gene, ".");print gene[3]":"$4"\t"$7}'"""\
              % (bed, bedgraph)
    sp.call(cmd, shell=True)

def main(counts, sites, sizes, bases):
    # tmpsizes = get_sizes(species)
    tmp = add_slop(sites, sizes, bases)
    map_counts(tmp, counts)
    os.remove(tmp)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("counts", help="read counts for a sample as bedgraph")
    p.add_argument("sites", help="merged polya sites as bed")
    p.add_argument("sizes", help="chromosome sizes for specific genome")
    # p.add_argument("-s", dest="species", default="hg18",
    #         help="ucsc table to use")
    p.add_argument("-b", dest="bases", type=int, default=2,
            help="increase region -b base pairs in each direction when \
            finding max read counts")
    args = vars(p.parse_args())
    main(**args)
