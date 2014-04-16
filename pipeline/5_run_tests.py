#!/usr/bin/env python
# encoding: utf-8
"""
Only the comparisons defined in the metadata sheet will be performed.

Run this from the command line -- Don't submit this to queue as is.
"""
import os
import sys
from bsub import bsub
from glob import glob
from string import Template
from toolshed import reader

results = "/vol1/home/brownj/projects/polya/results/common"
fisher_script = "/vol1/home/brownj/devel/polya/fisher_test.py"
fisher_results = "/vol1/home/brownj/projects/polya/results/common/fisher_results"
metadata = "/vol1/home/brownj/projects/polya/results/common/hub/metadata.tsv"
sample_column = "alias"
comparison_column = "comparisons"
files = glob("/vol1/home/brownj/projects/polya/results/common/*/*.counts.txt.gz")
fisher_submit = bsub("fisher", q="short", P="pillai_kabos_polya", verbose=True)

for t in reader(metadata, header=True):
    for compareto in t[comparison_column].split(","):
        for strand in ['pos', 'neg']:
            if len(compareto) == 0: continue
            a = t[sample_column]
            b = compareto
            fisher_result = "{fisher_results}/{a}_to_{b}.{strand}.fisher.txt.gz".format(fisher_results=fisher_results, a=a, b=b, strand=strand)

            if os.path.exists(fisher_result):
                print >>sys.stderr, ">> fisher comparison complete for", a, strand, "and", b, strand
                continue
            print >>sys.stderr, ">> fisher: comparing", a, "and", b, strand

            file_a = "{results}/{a}/{a}.{strand}.counts.txt.gz".format(results=results, a=a, strand=strand)
            file_b = "{results}/{b}/{b}.{strand}.counts.txt.gz".format(results=results, b=b, strand=strand)
            assert os.path.exists(file_a) and os.path.exists(file_b)

            fisher_cmd = ("python {fisher_script} {file_a} {file_b} | gzip -c > {fisher_result}").format(fisher_script=fisher_script, file_a=file_a, file_b=file_b, fisher_result=fisher_result)
            fisher_submit(fisher_cmd)
