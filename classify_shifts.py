#!/usr/bin/env python
# encoding: utf-8
"""
Characterize poly(A) shifts from DEXSeq test results as proximal (p) or
distal (d). Sample names are taken from the header (log2fold) of the DEXSeq
output files.
"""
import sys
import pandas as pd
from toolshed import reader
from itertools import groupby
from collections import OrderedDict

class StrandNotFound(Exception):
    pass

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

def gstrand(fname):
    if "pos" in fname:
        return "pos"
    if "neg" in fname:
        return "neg"
    raise StrandNotFound()

def grouper(iterable, col):
    for k, g in groupby(iterable, key=lambda t: t[col]):
        yield g

def shift(aid, afold, bid, bfold):
    """determine direction of change.
    pos to neg - proximal
    neg to pos - distal
    fold change does not flip - indicate whether both increased or decreased
    """
    if aid > bid:
        return shift(bid, bfold, aid, afold)
    else:
        assert afold != 0
        if afold < 0:
            if bfold > 0:
                return "distal"
            if bfold < 0:
                return "decreased"
        if afold > 0:
            if bfold < 0:
                return "proximal"
            if bfold > 0:
                return "increased"

def pairs(odict):
    elements = odict.items()
    for i in xrange(0, len(elements)-1):
        yield elements[i], elements[i+1]

def main(dexseq, pval, pval_cutoff):
    dex_runs = OrderedDict()
    for fname in dexseq:
        cols = reader(fname, header=False).next()[1:]
        try:
            a, b = sample_names(cols, "log2fold")
            strand = gstrand(fname)
        except StrandNotFound:
            print >>sys.stderr, ">> strand (pos, neg) must be in file names."
            sys.exit(1)
        except UnboundLocalError:
            print >>sys.stderr, ">> failed to get sample names for", fname
            print >>sys.stderr, ">> skipping..."
            continue
        log2fold = cols[-1]
        assert a != b
        run_id = "{a}_to_{b}.{strand}".format(**locals())
        dex_runs[run_id] = {}
        for group in grouper(reader(fname, header=True), "geneID"):
            results = OrderedDict()
            for site in group:
                try:
                    # p-value threshold filtering
                    if float(site[pval]) > pval_cutoff: continue
                except ValueError:
                    continue
                # fold change should be recorded from dexseq
                assert site[log2fold] != "NA"
                site_id = int(site['exonID'].rsplit(".")[-1])
                results[site_id] = {'fc':float(site[log2fold]), 'name':site['exonID'].lstrip('E')}
            if len(results) < 2: continue
            # iterating over the pairs involved in switching event
            for (aid, ad), (bid, bd) in pairs(results):
                # the direction of change
                direction = shift(aid, ad['fc'], bid, bd['fc'])
                comp = "{aname},{bname}".format(aname=ad['name'], bname=bd['name'])
                # complex name to ease creating multiindex dataframe
                dex_runs[run_id]["{gene}:{comp}".format(gene=site['geneID'], comp=comp)] = direction
    try:
        df = pd.DataFrame(dex_runs)
        # create multiindex via split
        df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], names=['Gene','Sites'])
        df.to_csv(sys.stdout, sep="\t", na_rep="na")
    except Exception:
        # empty dataframe
        print >>sys.stderr, "No significant sites were found."

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('dexseq', metavar="DEXSEQ", nargs="+",
            help="DEXSeq result files obtained via `run_dexseq.py`.")
    p.add_argument('-p', dest="pval_cutoff", default=0.05, type=float,
            help="p-value cutoff")
    p.add_argument('-v', dest="pval", default="padjust",
            choices=['pvalue, padjust'], help="which p used for cutoff")
    args = vars(p.parse_args())
    main(**args)
