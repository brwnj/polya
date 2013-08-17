#!/usr/bin/env python
# encoding: utf-8
"""
Fisher test per poly(a) site after normalizing by total count at the gene level.
Counts are divded by the total number of mapped reads to a particular gene and
multiplied by the mean total count across samples. q-value is only calculated on
genes that have been tested.
"""
import os
import sys
import rpy2
import bisect
import argparse
import tempfile
import numpy as np
import pandas as pd
import subprocess as sp
import pandas.rpy.common as com
from array import array
from itertools import combinations
from toolshed import reader, nopen
from rpy2.robjects.packages import importr

stats = importr('stats')

def qvality(pvals):
    """@brentp"""
    f = open(tempfile.mkstemp(suffix=".qvality")[1], "w")
    f.write("\n".join(map(str, pvals)))
    f.close()
    cmd = ['qvality', f.name]
    p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
    pvalues, peps, qvalues = [array('d') for _ in range(3)]
    for pmax, pep, q in (map(float, l.split("\t")) for i, l in enumerate(p.stdout) if i > 0):
        pvalues.append(pmax)
        peps.append(pep)
        qvalues.append(q)
    # make sure they are sorted.
    for pset in (pvalues, peps, qvalues):
        assert all(a >= b for a, b in zip(pset[1:], pset[:-1]))
    p.wait()
    if p.returncode != 255: # ?
        print >>sys.stderr, p.stderr.read(), p.returncode
    # delete temp file
    os.remove(f.name)
    return pvalues, peps, qvalues

def gqvalue(pvalue, pvalues, peps, qvalues):
    """@brentp"""
    # find the index among all p-values
    idx = bisect.bisect_left(pvalues, pvalue)
    # and return the posterior prob, q
    try:
        return pvalue, peps[idx], qvalues[idx]
    except IndexError:
        return pvalue, peps[idx - 1], qvalues[idx - 1]

def gsample(fname):
    return os.path.basename(fname).split(".")[0]

def gshift(lst):
    """
    increasing to decreasing    --->    proximal
    increasing to increasing    --->    no shift
    increasing to no change     --->    proximal

    decreasing to increasing    --->    distal
    decreasing to decreasing    --->    no shift
    decreasing to no change     --->    distal

    no change to decreasing     --->    proximal
    no change to increasing     --->    distal
    no change to no change      --->    no shift
    """
    res = []
    for (a, b) in lst:
        if a > b:
            res.append("decreasing")
        elif b > a:
            res.append("increasing")
        else:
            res.append("equal")
    assert len(res) == 2
    if res[0] == "equal":
        if res[1] == "decreasing":
            return "proximal"
        if res[1] == "increasing":
            return "distal"
        return "no shift"
    elif res[0] == "increasing":
        if res[1] == "decreasing":
            return "proximal"
        if res[1] == "increasing":
            return "no shift"
        return "proximal"
    else: # res[0] == "decreasing"
        if res[1] == "decreasing":
            return "no shift"
        if res[1] == "increasing":
            return "distal"
        return "distal"

def gfoldchange(df):
    # fix zero issues...
    if df.sum().sum() == 0:
        return 0
    df = df + 1
    df = df.div(df.min())
    # no sample to sample change was observed at either site
    if df.sum().sum() == 4:
        return 0
    return np.log2(df.where(df > 1).sum().sum() / len(df))

def _apply_qval(row, pvalues, peps, qvalues):
    row['q'] = gqvalue(row['p'], pvalues, peps, qvalues)[2]

def gcompression(fname):
    return "gzip" if fname.endswith(".gz") else None

def main(a, b):
    aid = gsample(a)
    bid = gsample(b)
    df = pd.read_table(a, header=None, names=["site", aid], index_col="site", compression=gcompression(a))
    tmp_df = pd.read_table(b, header=None, names=["site", bid], index_col="site", compression=gcompression(b))
    df = df.join(tmp_df)
    # create multiindex via split
    df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], names=['Gene','Site'])
    # unique genes
    genes = set([g for g in df.index.get_level_values('Gene')])
    res = {}    # fisher testing per site per gene
    pvals = []  # store the p-values for qvality
    for gene in genes:
        gs = df.ix[gene]
        # filter out all with only one site
        if len(gs) < 2: continue
        # flag genes without counts
        use_gene_for_q = gs.sum().any()
        # normalize and round down
        try:
            gs = (gs / gs.sum().astype('float') * gs.sum().mean()).astype('int')
        except:
            # don't normalize when one sample is all 0s
            pass
        # test each site pair across the gene
        for (sitea, siteb) in combinations(gs.index, 2):
            res["{gene}:{sitea}:{siteb}".format(gene=gene, sitea=sitea, siteb=siteb)] = {}
            ss = gs.ix[[sitea, siteb]]
            # transpose
            sst = ss.T
            # convert slice to r::dataframe
            rs = com.convert_to_r_dataframe(sst)
            # fisher exact
            p = stats.fisher_test(rs)[0][0]
            # for some reason 1 is rounding to slightly greater than 1
            p = 1 if p > 1 else p
            shift = gshift(sst.values)
            fc = gfoldchange(sst.astype("float"))
            res["{gene}:{sitea}:{siteb}".format(**locals())]["shift"] = shift
            res["{gene}:{sitea}:{siteb}".format(**locals())]["foldchange"] = fc
            res["{gene}:{sitea}:{siteb}".format(**locals())]["p"] = p
            res["{gene}:{sitea}:{siteb}".format(**locals())]["q"] = 1.0
            if use_gene_for_q: pvals.append(p)
    # calculate qvalues
    pvalues, peps, qvalues = qvality(pvals)
    # convert fisher results into dataframe
    fisherdf = pd.DataFrame(res).T
    fisherdf = fisherdf.ix[:, ["shift", "foldchange", "p", "q"]]
    fisherdf.index = pd.MultiIndex.from_tuples([x.split(":") for x in fisherdf.index], names=['Gene','SiteA','SiteB'])
    # update qvalue column with qvality output
    fisherdf.apply(_apply_qval, axis=1, **{"pvalues":pvalues, "peps":peps, "qvalues":qvalues})
    fisherdf.to_csv(sys.stdout, sep="\t", float_format="%.8g")

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("counts_a", help="first sample counts file in dexseq compatible format")
    p.add_argument("counts_b", help="second sample counts file")
    args = p.parse_args()
    main(args.counts_a, args.counts_b)
