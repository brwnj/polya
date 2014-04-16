#!/usr/bin/env python
# encoding: utf-8
"""
filter the site counts removing entries where gene is accounted for with same
count and same peak classification, but has multiple refseq annotations.
"""
from toolshed import reader
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main(count_file):
    observed = defaultdict(set)
    for toks in reader(count_file, header=['gene', 'fullname', 'count']):
        if "|" in toks['fullname']:

            gene_symbol = toks['fullname'].split("|")[1]

            if not observed.has_key(gene_symbol):
                observed[gene_symbol].add(toks['count'])
            elif toks['count'] in observed[gene_symbol]:
                continue

        print "{gene}\t{name}\t{count}".format(gene=toks['gene'], name=toks['fullname'], count=toks['count'])

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('counts')
    args = p.parse_args()
    main(args.counts)
