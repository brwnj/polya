#!/usr/bin/env python
# encoding: utf-8
"""
Given counts from 2 samples, generate fake replicates and run DEXSeq.
"""
import os
import sys
import tempfile
import subprocess
from bsub import bsub
from random import randint
from toolshed import reader
from itertools import combinations

class StrandMismatch(Exception):
    pass

class ComparisonComplete(Exception):
    pass

def replicate(fname, n=5):
    tmp = open(tempfile.mkstemp(suffix=".txt", dir=".")[1], 'w')
    for c in reader(fname, header=['name','count']):
        count = int(c['count'])
        low = count - n if count - n > 0 else 0
        count = randint(low, count + n)
        tmp.write("{name}\t{count}\n".format(name=c['name'], count=count))
    tmp.close()
    return tmp.name

def sample_name(fname):
    return os.path.basename(fname).split(".", 1)[0].split("_", 1)[0]

def get_strand(flst):
    strand = set(["neg" if "neg" in f else "pos" for f in flst])
    if not len(strand) == 1:
        raise StrandMismatch()
    return strand.pop()

def main(files, script, projid, queue, verbose):
    submit = bsub("dexseq", P=projid, q=queue, n="4", R="span[hosts=1]", verbose=verbose)
    for (a, b) in combinations(files, 2):
        result = ""
        try:
            strand = get_strand([a, b])
            sample_a = sample_name(a)
            sample_b = sample_name(b)
            result = "{sample_a}_vs_{sample_b}.{strand}.txt".format(**locals())
            out_file_check = os.path.splitext(result)[0]
            # this script is run in the dexseq_results directory
            for f in os.listdir("."):
                if f.startswith(os.path.splitext(result)[0]):
                    raise ComparisonComplete()
            if verbose:
                print >>sys.stderr, ">> comparing", sample_a, "and", sample_b
            rep_a = replicate(a)
            rep_b = replicate(b)
            cmd = ("Rscript {script} {sample_a},{sample_a}x {a},{rep_a} "
                    "{sample_b},{sample_b}x {b},{rep_b} "
                    "{result}").format(**locals())
            submit(cmd, job_cap=20)
        except StrandMismatch:
            continue
        except ComparisonComplete:
            print >>sys.stderr, ">> comparison complete:", result
            continue

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("script", help="full path to run_dexseq.R")
    p.add_argument("files", nargs="+", help="count files to be tested")
    req = p.add_argument_group("required arguments")
    req.add_argument("-p", dest="projid", required=True,
            help="project id for cluster usage tracking")
    p.add_argument("-q", dest="queue", default="normal",
            help="lsf queue")
    p.add_argument("-v", dest="verbose", action="store_true")
    args = vars(p.parse_args())
    main(**args)
