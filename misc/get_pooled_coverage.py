#!/usr/bin/env python
# encoding: utf-8
"""
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from bsub import bsub
from toolshed import reader
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def get_sample_name(fname, pattern):
    """
    >>>get_sample_name("PK96.pos.bedgraph.gz", ".bedgraph")
    'PK96.pos'
    """
    return fname.split(pattern)[0]

def get_compression_setting(fname):
    return "gzip" if fname.endswith(".gz") else None

def combinations(l):
    for i in xrange(1, len(l)):
        yield l[0], l[i]

def _sf_deseq(counts):
    """
    Calculate DESeq scaling factor per sample.
    """
    # masked array to discard inf, -inf, and nan
    ma = np.ma.masked_invalid(counts)
    return np.exp(np.ma.median(ma))

def norm_deseq(df):
    """
    Normalize by DESeq scaling factor, which is computed as the median of
    the ratio, for each row (gene), of its read count over its geometric
    mean across all samples. Return new counts dataframe.

    Details:
    --------
    http://genomebiology.com/2010/11/10/R106

    Parameters:
    -----------
    - df: pandas dataframe.
    """
    # log of counts
    lg = df.apply(np.log)
    # per sample: exponential(median(log(counts) - geometric mean))
    sf = lg.sub(lg.mean(axis=1), axis=0).apply(_sf_deseq, axis=0)
    # apply scaling
    df = df.div(sf, axis=1)
    return df

def norm_tc(df, verbose=False):
    """Normalize by total count."""
    # sum of all the counts
    sum_by_sample = df.sum()
    total_sum = sum_by_sample.sum()
    mean_total_count = float(total_sum / len(df.columns))
    if verbose:
        print >>sys.stderr, "Total Sum:", total_sum
        print >>sys.stderr, "Mean Total Count:", mean_total_count
        print >>sys.stderr, "Total Count by Sample:"
        sum_by_sample.to_csv(sys.stderr, sep="\t")
    # normalize each column
    df = df.apply(lambda x: (x / x.sum()) * mean_total_count)
    return df

def main(bedgraphs, metadata, verbose):
    pools = defaultdict(list)
    for toks in reader(metadata):
        for k, v in toks.iteritems():
            if k.startswith("Pool") and v == "TRUE":
                # get the samples
                pool_name = k.split("_")[-1]
                pools[pool_name].append(toks['alias'])

    for pool, samples in pools.iteritems():
        for strand in ["pos", "neg"]:

            if verbose:
                print >>sys.stderr, ">> processing", pool, strand

            files = [f for f in bedgraphs if os.path.basename(f).split(".")[0] in samples and strand in os.path.basename(f) and "UMI" not in f]
            if len(files) == 0: continue
            assert len(files) == len(samples), "All count files not found for {pool}".format(pool=pool)

            df_list = [pd.read_table(f, names=["chrom", "start", "stop", get_sample_name(f, ".bedgraph")], compression="gzip") for f in files]

            # fixing bedgraph files with spans greater than 1
            for df in df_list:
                df.stop = df.start + 1
                df.set_index(['chrom','start','stop'], inplace=True)

            # combine based on index
            combined_df = pd.concat(df_list, axis=1)
            # normalize the counts
            # normed_df = norm_deseq(combined_df)
            normed_df = norm_tc(combined_df, verbose)
            combined_df = None

            # round the normalized counts up to int
            normed_df = normed_df.apply(np.ceil)
            normed_df.fillna(0, inplace=True)
            # sum the rows
            normed_df[pool] = normed_df.sum(axis=1)
            # print results
            normed_df[pool].astype('int').to_csv("{pool}.{strand}.bedgraph".format(pool=pool, strand=strand), sep="\t")

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata")
    p.add_argument("bedgraphs", nargs="+", help="sorted sample bedgraph files")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()
    main(args.bedgraphs, args.metadata, args.verbose)
