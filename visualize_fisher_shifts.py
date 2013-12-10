#!/usr/bin/env python
# encoding: utf-8
"""
From the shifts and merged sites, create a bed12 for each comparison with a
trackline to visualize shift in UCSC browser. Unsorted output to stdout.
"""
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

class Bed(object):
    def __init__(self, toks):
        self.chrom = toks[0]
        self.start, self.stop = map(int, toks[1:3])
        self.name, self.score, self.strand = toks[3:6]

def bed12line(chrom, start, stop, strand, shift):
    score = 0
    color = "255,0,0"
    diff = stop - start - 1
    if strand == "+":
        arrow = "-" if shift == "proximal" else "+"
    if strand == "-":
        arrow = "+" if shift == "proximal" else "-"
    return [chrom, start, stop, shift, score, arrow, start, stop, color, 2, "1,1", "0,{diff}".format(diff=diff)]

def sites_to_dict(bed):
    d = {}
    for l in reader(bed, header=Bed):
        d[l.name] = l
    return d

def main(shifts, sites, cutoff=0.05, min_length=5):
    refsites = sites_to_dict(sites)
    for l in reader(shifts):
        if not l['shift'] == "proximal" and not l['shift'] == "distal": continue
        if float(l['q']) >= cutoff: continue
        sites = [l['SiteA'], l['SiteB']]
        # sort by ascending site IDs
        sites.sort(key=lambda x: int(x.split(".")[-1]))
        downstream, upstream = sites

        a = refsites[downstream]
        b = refsites[upstream]

        if abs(a.start - b.start) < min_length: continue
        # if the strand is negative, sites are classified in descending order
        if a.strand == "-":
            print "\t".join(map(str, bed12line(a.chrom, b.start, a.stop, a.strand, l['shift'])))
        else:
            print "\t".join(map(str, bed12line(a.chrom, a.start, b.stop, a.strand, l['shift'])))

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("shifts", help="output of `fisher_test.py`")
    p.add_argument("sites", help="output of `merge_sites.py`")
    p.add_argument("-q", dest="cutoff", default=0.05, type=float, help="suppress peaks with a q-value less than cutoff")
    p.add_argument("-l", dest="min_length", default=5, type=int, help="suppress peaks with a distance between them less than min_length")
    args = p.parse_args()
    main(args.shifts, args.sites, args.cutoff, args.min_length)
