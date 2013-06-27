#!/usr/bin/env python
# encoding: utf-8
"""
Given counts from 2 samples, generate fake replicates and run DEXSeq.
"""
import os, tempfile, subprocess
from bsub import bsub
from random import randint
from toolshed import reader
from itertools import combinations

class StrandMismatch(Exception):
    pass

def replicate(fname, n=5):
    tmp = open(tempfile.mkstemp(suffix=".txt")[1], 'w')
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

def main(files, script, projid, queue):
    submit = bsub("dexseq", P=projid)
    jobids = []
    for (a, b) in combinations(files, 2):
        try:
            strand = get_strand([a, b])
            sample_a = sample_name(a)
            sample_b = sample_name(b)
            rep_a = replicate(a)
            rep_b = replicate(b)
            cmd = ("Rscript {script} {sample_a},{sample_a}x {a},{rep_a} "
                    "{sample_b},{sample_b}x {b},{rep_b} "
                    "{sample_a}_vs_{sample_b}.{strand}.txt").format(**locals())
            wait = submit(cmd)
            cmd = "rm {rep_a} {rep_b}".format(**locals())
            jobids.append(bsub("dexseq_cleanup", P=projid, w=wait)(cmd))
        except StrandMismatch:
            continue
    # bsub.poll(jobids)
    # check the output files
    # resubmit anything with all NA

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
    args = vars(p.parse_args())
    main(**args)
