#!/usr/bin/env python
# encoding: utf-8
"""
"""

def bed12line(chrm, strt, stp, shft):
    arrw = "+" if shft == ""
    l = int(strt) - int(stp)
    return "{chrm}\t{strt}\t{stp}\t{shft}\t0\t{arrw}\t{strt}\t{stp}\t0\t1\t{l}\t0".format(**locals())

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
