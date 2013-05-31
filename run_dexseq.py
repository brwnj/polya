#!/usr/bin/env python
# encoding: utf-8
"""
Given counts from 2 samples, generate a fake replicate and run DEXSeq.
"""
import os, tempfile, subprocess
from random import randint
from toolshed import reader
from itertools import combinations

class StrandMismatch(Exception):
    pass

def replicate(fname, n=5):
    tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
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

def dexseq(a, b, script):
    strand = get_strand([a, b])
    sample_a = sample_name(a)
    sample_b = sample_name(b)
    rep_a = replicate(a)
    rep_b = replicate(b)
    cmd = ("Rscript {script} {sample_a},{sample_a}x {a},{rep_a} "
            "{sample_b},{sample_b}x {b},{rep_b} "
            "{sample_a}_vs_{sample_b}.{strand}.txt").format(**locals())
    subprocess.call(cmd, shell=True)
    os.remove(rep_a)
    os.remove(rep_b)

def main(files, script):
    for (a, b) in combinations(files, 2):
        try:
            dexseq(a, b, script)
        except StrandMismatch:
            continue

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("script", help="run_dexseq.R")
    p.add_argument("files", nargs="+", help="count files to be tested")
    args = vars(p.parse_args())
    main(**args)
