#!/usr/bin/env python
# encoding: utf-8
"""
Fisher test per poly(a) site after normalizing by total count at the gene level.
Counts are divded by the total number of mapped reads to a particular gene and
multiplied by the mean total count across samples.

q-value is only calculated on genes that have been successfully tested.
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

def get_qvalue(pvalue, pvalues, peps, qvalues):
    """@brentp"""
    # find the index among all p-values
    idx = bisect.bisect_left(pvalues, pvalue)
    # and return the posterior prob, q
    try:
        return pvalue, peps[idx], qvalues[idx]
    except IndexError:
        return pvalue, peps[idx - 1], qvalues[idx - 1]

def main(a, b):
    # map counts to polya sites
    # count_files = get_count_files(args.bed, args.bedgraphs)
    # build the dataframe
    
    # Total count (TC): Gene counts are divided by the total number of mapped
    # reads (or library size) associated with their lane and multiplied by the
    # mean total count across all the samples of the dataset.
    
    # will have to split the first column into something meaningful
    df = pd.read_csv...(a...)
    # read_table(cfile, header=None, index_col=[0,1], names=["gene", "polya", sample])
    temp_df = pd.read_csv...(b...)
    # maybe split into multiindex at this point
    df = df.join(temp_df)

    # unique genes
    genes = set()
    for gene in count_df.index.get_level_values('gene'):
        genes.add(gene)

    # fisher testing per site per gene
    fisher_test = {'p_value':{}}
    # store the p-values for qvality
    pvals = []
    for i, gene in enumerate(genes, start=1):
        gene_slice = count_df.ix[gene]

        # more than one polya site and not all zero
        if len(gene_slice) < 2 or gene_slice.sum().any() == 0: continue

        # convert slice to r::dataframe
        r_slice = com.convert_to_r_dataframe(gene_slice.transpose())

        # fisher exact
        p = stats.fisher_test(r_slice)[0][0]

        # for some reason 1 is rounding to slightly greater than 1
        if p > 1: p = 1
        fisher_test["p_value"][gene] = p
        pvals.append(p)

    # double check
    assert all(0 <= p <= 1 for p in pvals)

    # calculate qvalues
    pvalues, peps, qvalues = qvality(pvals)
    # convert dict of dicts to dataframe
    fisher_df = pd.DataFrame(fisher_test)
    # remove multiindex and set to gene level
    df = df.reset_index().set_index("gene")
    # add p-value column
    df = df.join(fisher_df)
    df = df.reset_index()
    # print everything including pvals to temp
    pvals_f = tempfile.mkstemp(suffix=".pvals")[1]
    df.to_csv(pvals_f, sep="\t", header=False, index=False, na_rep="1.0", float_format="%.6g")
    # append pep and qvalue to output
    for line in nopen(pvals_f):
        toks = line.rstrip("\r\n").split("\t")
        p, pep, q = get_qvalue(float(toks[4]), pvalues, peps, qvalues)
        print "%s\t%.6g\t%.6g" % (line.rstrip("\r\n"), pep, q)
    os.remove(pvals_f)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("counts_a", help="first sample counts file in dexseq compatible format")
    p.add_argument("counts_b", help="second sample counts file")
    args = p.parse_args()
    main(args.counts_a, args.counts_b)
