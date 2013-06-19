#!/usr/bin/env python
# encoding: utf-8
"""
From the shifts and merged sites, create a bed12 in order to improve shift
visualization in UCSC browser.
"""
from toolshed import reader

class Bed(object):
    def __init__(self, toks):
        self.chrom = toks[0]
        self.start, self.stop = map(int, toks[1:3])
        self.name, self.score, self.strand = toks[3:6]

def bed12line(chr, start, stop, strand, shift):
    if strand == "+":
        arrow = "-" if shift == "proximal" else "+"
    if strand == "-":
        arrow = "+" if shift == "proximal" else "-"
    score = 0
    color = "255,0,0"
    diff = stop - start - 1
    return ("{chr}\t{start}\t{stop}\t{shift}\t{score}\t{arrow}\t"
            "{start}\t{stop}\t{color}\t2\t1,1\t0,{diff}").format(**locals())

def sites_to_dict(bed):
    d = {}
    for l in reader(bed, header=Bed):
        d[l.name] = l
    return d

def main(args):
    sites = sites_to_dict(args.sites)
    cols = reader(args.shifts, header=False).next()
    shift_col = cols[-1]
    for l in reader(args.shifts):
        a, b = l['Sites'].split(",")
        a = sites[a]
        b = sites[b]
        shift = l[shift_col]
        if shift != "proximal" and shift != "distal": continue
        print bed12line(a.chrom, a.start, b.stop, a.strand, shift)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("shifts", help="output of `classify_shifts.py`")
    p.add_argument("sites", help="output of `merge_sites.py`")
    args = p.parse_args()
    main(args)
