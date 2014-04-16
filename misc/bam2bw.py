#!/usr/bin/env python
# encoding: utf-8
"""
Convert BAM to BW. If BW is not found in the directory of the bam, this script
will overwrite existing intermediate files in its creation.
"""

import os
import sys
import pysam
from bsub import bsub
from contextlib import contextmanager
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


@contextmanager
def indexed_bam(bam):
    assert os.path.exists(bam)
    if not os.path.exists(bam + ".bai"):
        print >>sys.stderr, "indexing", bam
        pysam.index(bam)
    reader = pysam.Samfile(bam)
    try:
        yield reader
    finally:
        reader.close()


def main(bam, project_id, fprime, queue, verbose):
    """
    converts /vol/home/sample.bam to:

        /vol/home/user/sample_pos.bedgraph.gz
        /vol/home/user/sample_neg.bedgraph.gz
        /vol/home/user/sample_pos.bw
        /vol/home/user/sample_neg.bw

    """

    sample, ext = bam.split(".bam", 1)
    sizes_file = sample + ".sizes"

    with indexed_bam(bam) as bamfh, open(sizes_file, 'w') as sizesfh:
        for chrom, length in zip(bamfh.references, bamfh.lengths):
            print >>sizesfh, "{chrom}\t{length}".format(**locals())

    submit = bsub("bam2bg", P=project_id, q=queue, verbose=verbose)

    # convert bam to stranded bg then bw
    for symbol, strand in zip(["+", "-"], ["pos", "neg"]):

        # check for existing bw
        bigwig = "{sample}_{strand}.bw".format(**locals())
        if exists(bigwig) and getsize(bigwig) > 0: continue

        bedgraph = "{sample}_{strand}.bedgraph".format(**locals())

        bam_to_bg = ["bedtools", "genomecov", "-strand", symbol, "-bg", "-ibam", bam]
        if fprime:
            bam_to_bg.append("-5")
        bam_to_bg.append("| bedtools sort -i - > {bedgraph}".format(**locals()))

        bg_to_bw = "bedGraphToBigWig {bedgraph} {sizes} {bigwig}".format(**locals())
        gzip_bg = "gzip -f {bedgraph}".format(**locals())

        submit(" ".join(bam_to_bg)).then(bg_to_bw, "bg2bw").then(gzip_bg, "gzipbg")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("bam", metavar="BAM", help="alignment file")
    p.add_argument("-5", dest="fprime", action="store_true", help="5' read counts")
    p.add_argument("-p", dest="project_id", default="bam2bw", help="can specify for usage tracking on cluster")
    p.add_argument("-q", dest="queue", default="normal", help="lsf queue")
    p.add_argument("-v", dest="verbose", action="store_true", help="print job commands")
    args = vars(p.parse_args())
    main(**args)
