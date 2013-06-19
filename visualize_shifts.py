#!/usr/bin/env python
# encoding: utf-8
"""
"""

def bed12line(chr, start, stop, shift):
    # needs to account for strand of gene and shift direction
    arrow = "+"
    score = 0
    color = "255,0,0"
    diff = int(stop) - int(start) - 1
    return "{chr}\t{start}\t{stop}\t{shift}\t{score}\t{arrow}\t{start}\t{stop}\t{color}\t2\t1,1\t0,{diff}".format(**locals())

def main(args):
    

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    shifts - output of `classify_shifts.py`
    merged sites - output of `merge_sites.py`
    p.add_argument()
    args = p.parse_args()
    main(args)
