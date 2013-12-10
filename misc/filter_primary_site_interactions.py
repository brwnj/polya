#!/usr/bin/env python
# encoding: utf-8
"""
Filter Fisher test results to show only those involving a primary site from
either sample in the comparison. Expects file names coming in to be
sample1_to_sample2.strand.etc.txt.
"""
import os
import pandas
import signal
from toolshed import header, reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def get_ab(fname):
    """
    >>> get_ab("NBT_to_TS-ERP.neg.fisher.txt.gz")
    ('NBT.neg', 'TS-ERP.neg')
    >>> get_ab("PK95_to_PK105.pos.fisher.txt.gz")
    ('PK95.pos', 'PK105.pos')
    """
    fname = os.path.basename(fname)
    parts = fname.split(".", 1)[0].split("_")
    a = parts[0]
    b = parts[2]
    strand = "pos" if "pos" in fname else "neg"
    return a + "." + strand, b + "." + strand

def get_comp(fname):
    return "gzip" if fname.endswith(".gz") else None

def main(tests, xref, head):
    df = pandas.read_table(xref, index_col=[0], compression=get_comp(xref))
    a, b = get_ab(tests)
    out_order = header(tests)
    if head:
        print "\t".join(out_order)
    for l in reader(tests):
        try:
            if df[a][l['SiteA']] or df[b][l['SiteA']] or df[a][l['SiteB']] or df[b][l['SiteB']]:
                print "\t".join(l[h] for h in out_order)
        except KeyError:
            # definitely not a primary site
            pass

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('test_results', help="fisher test result text file")
    p.add_argument('primary_site_table', help="primary site annotation table")
    p.add_argument('--header', dest="head", action='store_true', help="maintain the file header")
    args = p.parse_args()
    main(args.test_results, args.primary_site_table, args.head)
