#!/usr/bin/env python
# encoding: utf-8
"""
From the shifts and merged sites, create a bed12 for each comparison with a
trackline to visualize shift in UCSC browser. Strand is determined from the file
name and must be present.
"""
import sys
import operator
from toolshed import reader

class StrandNotFound(Exception):
    pass

class Bed(object):
    def __init__(self, toks):
        self.chrom = toks[0]
        self.start, self.stop = map(int, toks[1:3])
        self.name, self.score, self.strand = toks[3:6]

def get_strand(fname):
    if "pos" in fname:
        return "pos"
    if "neg" in fname:
        return "neg"
    raise StrandNotFound()

def bed12line(chrom, start, stop, strand, shift):
    if strand == "+":
        arrow = "-" if shift == "proximal" else "+"
    if strand == "-":
        arrow = "+" if shift == "proximal" else "-"
    score = 0
    color = "255,0,0"
    diff = stop - start - 1
    return [chrom, start, stop, shift, score, arrow,
                start, stop, color, 2, "1,1", 0, diff]

def sites_to_dict(bed):
    d = {}
    for l in reader(bed, header=Bed):
        d[l.name] = l
    return d

def shifts_to_dict(cols, fname):
    d = {}
    for c in cols:
        d[c] = {}
        for l in reader(fname):
            if not l[c] == "proximal" and not l[c] == "distal": continue
            d[c][l['Sites']] = l[c]
    return d

def main(shifts, sites):
    try:
        strand = get_strand(shifts)
    except StrandNotFound:
        print >>sys.stderr, "\nStrand ('pos', 'neg') must be in file name.\n"
        sys.exit(1)
    refsites = sites_to_dict(sites)
    cols = reader(shifts, header=False).next()
    comparisons = cols[2:]
    shifts_d = shifts_to_dict(comparisons, shifts)
    for comparison, all_sites in shifts_d.iteritems():
        lines = []
        for (site, shift) in all_sites.iteritems():
            a, b = site.split(",")
            a = refsites[a]
            b = refsites[b]
            lines.append(bed12line(a.chrom, a.start, b.stop, a.strand, shift))
        lines = sorted(lines, key=operator.itemgetter(0, 1))
        f = open("{comparison}.{strand}.bed".format(**locals()), 'wb')
        for line in lines:
            print >>f, "\t".join(map(str, line))
        f.close()

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("shifts", help="output of `classify_shifts.py`")
    p.add_argument("sites", help="output of `merge_sites.py`")
    args = p.parse_args()
    main(args.shifts, args.sites)
