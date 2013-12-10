#!/usr/bin/env python
# encoding: utf-8
"""
Since counts are in the form of bedgraph, counts as a result of this script
are unstranded. Any input should be strand specific.
"""
from toolshed import nopen
from operator import attrgetter
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

class Bed(object):
    __slots__ = ['chrom','start','stop','name','score','strand','mapped_score']
    def __init__(self, args):
        # clean up the bedtools input
        args = args.strip().split("\t")
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        self.start = int(self.start)
        self.stop = int(self.stop)
        self.score = int(self.score)
        self.mapped_score = int(self.mapped_score)

    def __repr__(self):
        return "BED({chr}:{start}-{stop})".format(chr=self.chrom, start=self.start, stop=self.stop)

    @property
    def gene(self):
        return self.name.split(".", 1)[1].split("|")[0].split(".")[0]

    @property
    def polya(self):
        return self.name.rsplit(".", 1)[1]

def main(counts, sites):
    cmd = "|bedtools map -c 4 -o max -null 0 -a {bed} -b {bedgraph}".format(bed=sites, bedgraph=counts)
    lines = [Bed(line) for line in nopen(cmd) if line.strip()]
    for bed in sorted(lines, key=attrgetter('gene', 'polya')):
        print "{gene}\t{name}\t{count}".format(gene=bed.gene, name=bed.name, count=bed.mapped_score)

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("counts", help="stranded counts as bedgraph")
    p.add_argument("sites", help="stranded polya sites as bed. most likely you'll have added slop onto original entries.")
    args = p.parse_args()
    main(args.counts, args.sites)
