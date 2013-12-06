#!/usr/bin/env python
# encoding: utf-8
"""
Fisher test per poly(a) site after normalizing by total count at the gene level.
Counts are divded by the total number of mapped reads to a particular gene and
multiplied by the mean total count across samples. q-value is only calculated on
genes that have been tested. Input file format is Gene, Polya Site Name, and
Count separated by tabs.
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
from rpy2.robjects.packages import importr
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

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

def get_qvalue(pvalue, pvalues, peps, qvalues):
    """@brentp"""
    # find the index among all p-values
    idx = bisect.bisect_left(pvalues, pvalue)
    # and return the posterior prob, q
    try:
        return pvalue, peps[idx], qvalues[idx]
    except IndexError:
        return pvalue, peps[idx - 1], qvalues[idx - 1]

def get_sample_name(fname):
    return os.path.basename(fname).split(".")[0]

def get_shift_direction(count_matrix):
    """
    slope...
    
    increasing to decreasing    --->    proximal
    increasing to increasing    --->    no shift
    increasing to no change     --->    proximal

    decreasing to increasing    --->    distal
    decreasing to decreasing    --->    no shift
    decreasing to no change     --->    distal

    no change to decreasing     --->    proximal
    no change to increasing     --->    distal
    no change to no change      --->    no shift
    
    >>> import numpy as np
    >>> arr = np.ndarray(shape=(2,2), dtype=int, buffer=np.array(0,2,16,3))
    >>> get_shift_direction(arr)
    'proximal'
    """
    assert len(count_matrix) == 2
    shift_types = {'equal':{'decreasing':'proximal',
                            'increasing':'distal',
                            'equal':'no shift'},
                   'increasing':{'decreasing':'proximal',
                                 'increasing':'no shift',
                                 'equal':'proximal'},
                   'decreasing':{'decreasing':'no shift',
                                  'increasing':'distal',
                                  'equal':'distal'}}
    observations = []
    for (a, b) in count_matrix:
        if a > b:
            observations.append("decreasing")
        elif b > a:
            observations.append("increasing")
        else:
            observations.append("equal")
    return shift_types[observations[0]][observations[1]]

def get_fold_change(df):
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
    row['q'] = get_qvalue(row['p'], pvalues, peps, qvalues)[2]

def get_compression_setting(fname):
    return "gzip" if fname.endswith(".gz") else None

def get_dataframe(count_file, sample_id):
    comp = get_compression_setting(count_file)
    df = pd.read_table(count_file, header=None, names=["gene", "site", sample_id],
                        index_col=["gene", "site"], compression=comp)
    return df

def main(a, b):
    aid = get_sample_name(a)
    bid = get_sample_name(b)
    df = get_dataframe(a, aid)
    tmp_df = get_dataframe(b, bid)
    df = df.join(tmp_df)
    # unique genes
    genes = set([g for g in df.index.get_level_values('gene')])
    # fisher testing per site per gene
    res = {}
    # store the p-values for qvality
    pvals = []
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
        # order the index based on site position
        ordered_index = sorted(gs.index, key=lambda x: int(x.rsplit(".", 1)[-1]))
        # test each site pair across the gene
        for (sitea, siteb) in combinations(ordered_index, 2):
            # the sites must be in order for the shift direction to be correct
            assert int(sitea.rsplit(".", 1)[-1]) < int(siteb.rsplit(".", 1)[-1])
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
            shift = get_shift_direction(sst.values)
            fc = get_fold_change(sst.astype("float"))
            res["{gene}:{sitea}:{siteb}".format(gene=gene, sitea=sitea, siteb=siteb)]["shift"] = shift
            res["{gene}:{sitea}:{siteb}".format(gene=gene, sitea=sitea, siteb=siteb)]["foldchange"] = fc
            res["{gene}:{sitea}:{siteb}".format(gene=gene, sitea=sitea, siteb=siteb)]["p"] = p
            res["{gene}:{sitea}:{siteb}".format(gene=gene, sitea=sitea, siteb=siteb)]["q"] = 1.0
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
    p = ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("counts_a", help="count file as described in docstring")
    p.add_argument("counts_b", help="second count file with same format")
    args = p.parse_args()
    main(args.counts_a, args.counts_b)
