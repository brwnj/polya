#!/usr/bin/env python
# encoding: utf-8
"""
Convert bed to bigbed. Score field will be scaled to 1000 as max. UCSC
bedToBigBed >= v. 2.5 must be in PATH.
"""
import os
import sys
import argparse
import tempfile
import subprocess as sp
from toolshed import reader


class Bed(object):
    __slots__ = ['chrom','start','stop','name', 'score', 'strand', 'thickstart',
            'thickstop', 'itemrgb', 'blockcount', 'blocksizes', 'blockstarts']

    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        self.score = int(round(float(self.score)))

    def __str__(self):
        fields = []
        for s in self.__slots__:
            try:
                fields.append(getattr(self, s))
            except AttributeError:
                continue
        return "\t".join(fields)


def scale(score, omax, omin, smax, smin):
    """
    >>> scale(2871, 4871, 0, 1000, 0)
    589
    """
    try:
        return ((smax - smin) * (score - omin) / (omax - omin)) + smin
    except ZeroDivisionError:
        return 0


def main(bed, sizes, mx, mn, bedtype):
    if os.path.getsize(bed) < 1:
        print >>sys.stderr, ">> skipping empty file", bed
        sys.exit(0)
    name = bed.split(".bed")[0]
    scores = []

    # find max and min
    for b in reader(bed, header=Bed):
        scores.append(b.score)
    try:
        scoremx = float(max(scores))
    except ValueError:
        print >>sys.stderr, bed, "appears empty."
        sys.exit(1)
    scoremn = float(min(scores))
    scaledtmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')

    # scale the values
    for b in reader(bed, header=Bed):
        b.score = str(int(scale(b.score, scoremx, scoremn, mx, mn)))
        print >>scaledtmp, b
    scaledtmp.close()

    # make sure the bed is sorted
    sortedtmp = tempfile.mkstemp(suffix=".bed")[1]
    cmd = "sort -k1,1 -k2,2n {scaledtmp} > {sortedtmp}".format(
            scaledtmp=scaledtmp.name, sortedtmp=sortedtmp)
    sp.call(cmd, shell=True)

    # convert to bb
    cmd = "bedToBigBed -type={bedtype} {tmp} {sizes} {name}.bb".format(
            tmp=sortedtmp, sizes=sizes, name=name, bedtype=bedtype)
    sp.call(cmd, shell=True)
    os.remove(scaledtmp.name)
    os.remove(sortedtmp)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("sizes", metavar="SIZES", help="genome sizes file")
    p.add_argument("bed", metavar="BED",
            help="unzipped bed6 to convert to bigbed")
    p.add_argument("--score-max", dest="scoremx", default=1000, type=int,
            help="maximum score to be observed")
    p.add_argument("--score-min", dest="scoremn", default=0, type=int,
            help="mininum score to be observed")
    p.add_argument("--type", dest="bedtype", default="bed6",
            help="examples: -type=bed6 or -type=bed6+ or -type=bed6+3")
    args = p.parse_args()
    if 0 < args.scoremx > 1000 or 0 > args.scoremn > 1000 \
            or args.scoremn > args.scoremx:
        print >>sys.stderr, "Scores need to satisfy: 0 <= input_val <= 1000"
        sys.exit(1)
    main(args.bed, args.sizes, args.scoremx, args.scoremn, args.bedtype)
