#!/usr/bin/env python
# encoding: utf-8
"""
Characterize poly(a) shifts from DEXSeq test results. Sample names are taken
from the header (log2fold) of the DEXSeq output files.

d  - distal shift
nt - notest
ns - not significant
p  - proximal shift

"""
import sys
from itertools import groupby
from toolshed import reader

def parse_name(name):
    """
    >>> parse_name("MP51MP51x")
    'MP51'
    """
    return name[:len(name)/2]

def sample_names(cols, item):
    for c in cols:
        if not item in c: continue
        b, a = c.rstrip(")").split("(", 1)[-1].split("/")
    return parse_name(a), parse_name(b)

def grouper(iterable, col):
    for k, g in groupby(iterable, key=lambda t: t[col]):
        yield g

def shift(fld_changes):
    """determine direction of change."""    
    prev = ""
    direction = ""
    lst = [fld_changes[x] for x in sorted(fld_changes.keys())]
    for v in lst:
        if not prev:
            prev = "pos" if v > 0 else "neg"
            continue
        if v > 0 and prev == "pos": continue # still a proximal shift
        if v > 0 and prev == "neg":
            direction += "d"
            prev = "pos"
        if v < 0 and prev == "neg": continue # still a distal shift
        if v < 0 and prev == "pos":
            direction += "p"
            prev = "neg"
    return direction

def main(dexseq, pval):
    dex_runs = {}
    for fname in dexseq:
        cols = reader(fname, header=False).next()[1:]
        a, b = sample_names(cols, "log2fold")
        log2fold = cols[-1]
        assert a != b
        run_id = "[{a},{b}]".format(**locals())
        dex_runs[run_id] = {}
        for group in grouper(reader(fname, header=True), "geneID"):
            results = {}
            for site in group:
                site_num = int(site['exonID'].rsplit(".")[-1])
                # NA for any gene without multiple sites
                if site['padjust'] == "NA": continue
                # p-value threshold filtering
                if float(site['padjust']) > pval: continue
                results[site_num] = float(site[log2fold])
            if len(results) < 2: continue
            # determine direction of change in other cases
            dex_runs[run_id][site['geneID']] = shift(results)
    print dex_runs
    # into datatable with na vals
    # print to facility gene search in viewing
    # desired output is something like:
    # gene    MP55_to_MP56    MP56_to_MP57
    # DEK     p               p
    # MALAT1  d               d

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('dexseq', metavar="DEXSEQ", nargs="+",
            help="DEXSeq results files via `run_dexseq.py`.")
    p.add_argument('-p', dest="pval", default=0.05, type=float,
            help="p-value cutoff")
    args = vars(p.parse_args())
    main(**args)
