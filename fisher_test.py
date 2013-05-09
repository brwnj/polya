#!/usr/bin/env python
# encoding: utf-8
"""
fisher per gene over poly(a) sites. input is 2 or more bedgraphs to compare.

q-value is only calculated on genes that have been successfully tested.
"""
import os
import sys
import rpy2
import bisect
import tempfile
import numpy as np
import pandas as pd
import os.path as op
import subprocess as sp
import pandas.rpy.common as com
from array import array
from rpy2.robjects.packages import importr
from toolshed import reader, nopen

stats = importr('stats')

def get_count_files(bed, bedgraphs):
    samples = {}
    for bg in bedgraphs:
        sample_id = op.basename(bg).split(".")[0]
        print >> sys.stderr, ">> converting %s" % sample_id
        tf = open(tempfile.mkstemp(suffix=".count")[1], 'w')
        # map the counts onto the regions
        # only want the peak max here, not the sum over the slop
        cmd = "|bedtools map -c 4 -o max -null 0 -a %s -b %s" % (bed, bg)
        result_header = "chrom start stop name score strand count".split()
        # only print out gene, polya, and the count data
        for r in reader(cmd, header=result_header):
            gene, polya = r['name'].split(":")
            tf.write("%s\t%s\t%s\n" % (gene, polya, r['count']))
        tf.close()
        samples[sample_id] = tf.name
    return samples

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

def main(args):
    if len(args.bedgraphs) < 2:
        print >>sys.stderr, ">> at least 2 bedgraphs..."
        sys.exit(1)
    # map counts to polya sites
    count_files = get_count_files(args.bed, args.bedgraphs)
    # build the dataframe
    count_df = pd.DataFrame()
    for sample, cfile in count_files.iteritems():
        print >> sys.stderr, ">> reading %s..." % sample
        if count_df.empty:
            count_df = pd.io.parsers.read_table(cfile, header=None, \
                        index_col=[0,1], names=["gene", "polya", sample])
        else:
            temp = pd.io.parsers.read_table(cfile, header=None, \
                    index_col=[0,1], names=["gene", "polya", sample])
            # append the column onto table
            count_df = count_df.join(temp)
        os.remove(cfile)
    # unique genes
    genes = set()
    for gene in count_df.index.get_level_values('gene'):
        genes.add(gene)
    # fisher testing over genes
    print >> sys.stderr, ">> processing..."
    fisher_test = {}
    fisher_test['p_value'] = {}
    # store the p-values for qvality
    pvals = []
    for i, gene in enumerate(genes, start=1):
        if i % 5000 == 0:
            print >> sys.stderr, ">> processed %d genes..." % i
        gene_slice = count_df.ix[gene]
        # more than one polya site and not all zero
        if len(gene_slice) < 2 or gene_slice.sum().any() == 0: continue
        # convert slice to r::dataframe
        r_slice = com.convert_to_r_dataframe(gene_slice.transpose())
        # fisher exact
        # try
        p = stats.fisher_test(r_slice, workspace=20000000)[0][0]
        # for some reason 1 is rounding to slightly greater than 1
        # except: increment workspace
        if p > 1: p = 1
        fisher_test["p_value"][gene] = p
        pvals.append(p)
    # double check
    assert all(0 <= p <= 1 for p in pvals)
    # run qvality
    pvalues, peps, qvalues = qvality(pvals)
    # convert dict of dicts to dataframe
    fisher_test_df = pd.DataFrame(fisher_test)
    # remove multiindex and set to gene level
    count_df = count_df.reset_index().set_index("gene")
    # add p-value column
    count_df = count_df.join(fisher_test_df)
    count_df = count_df.reset_index()
    # print everything including pvals to temp
    pvals_f = tempfile.mkstemp(suffix=".pvals")[1]
    count_df.to_csv(pvals_f, sep="\t", header=False, index=False, na_rep="1.0", float_format="%.6g")
    # append pep and qvalue to output
    for line in nopen(pvals_f):
        toks = line.rstrip("\r\n").split("\t")
        p, pep, q = get_qvalue(float(toks[4]), pvalues, peps, qvalues)
        print "%s\t%.6g\t%.6g" % (line.rstrip("\r\n"), pep, q)
    os.remove(pvals_f)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bed", help="reference bed of polya sites. the name field is gene:polya.")
    p.add_argument("bedgraphs", nargs="+", help="bedgraphs of samples to compare.")
    main(p.parse_args())