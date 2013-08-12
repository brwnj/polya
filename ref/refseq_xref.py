#!/usr/bin/env python
# encoding: utf-8
"""
rename exons to gene symbol using refseq xref and refseq bed of exons.
"""
import sys
import argparse
from toolshed import reader

class Bed(object):
    __slots__ = ['chrom','start','stop','name','score','strand']
    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)

    def __str__(self):
        return "\t".join([getattr(self, s) for s in self.__slots__])

def xref_to_dict(fname, a, b):
    d = {}
    for l in reader(fname):
        d[l[a]] = l[b]
    return d

def main(exons, xref):
    refseq_xref = xref_to_dict(xref, "From", "To")
    for b in reader(exons, header=Bed):
        let, num, _junk = b.name.split("_", 2)
        refseq_id = "{let}_{num}".format(**locals())
        try:
            b.name = refseq_xref[refseq_id]
        except KeyError:
            b.name = refseq_id
        print b

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('exons')
    p.add_argument('xref')
    args = p.parse_args()
    main(args.exons, args.xref)
