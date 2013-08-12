#!/usr/bin/env python
# encoding: utf-8
"""
exon reference to whole gene.
"""
import argparse
from toolshed import reader

class Bed(object):
    __slots__ = ['chrom','start','stop','name','score','strand']
    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        self.start = int(self.start)
        self.stop = int(self.stop)

    def __repr__(self):
        return "Bed({chr}:{name})".format(chr=self.chrom, name=self.name)

    def __str__(self):
        return "\t".join([str(getattr(self, s)) for s in self.__slots__])

def main(bed):
    mnd = {}
    mxd = {}
    for b in reader(bed, header=Bed):
        try:
            if mnd[b.name] > b.start:
                mnd[b.name] = b.start
            if mxd[b.name] < b.stop:
                mxd[b.name] = b.stop
        except KeyError:
            mnd[b.name] = b.start
            mxd[b.name] = b.stop
    seen = set()
    for b in reader(bed, header=Bed):
        if b.name in seen: continue
        seen.add(b.name)
        b.start = mnd[b.name]
        b.stop = mxd[b.name]
        print b

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('bed')
    args = p.parse_args()
    main(args.bed)
