#!/usr/bin/env python
# encoding: utf-8
"""
From sample counts, give table of Gene, Site, and 0/1 for a site being the 
primary site for any given gene.
"""
import os
import sys
import pandas
import argparse
from toolshed import reader
from itertools import groupby

class CountRecord(object):
    def __init__(self, args):
        self.gene_site = args[0]
        self.count = int(args[1])
    
    def __repr__(self):
        return "CountRecord:{site}".format(site=self.site)

    @property
    def gene(self):
        return self.gene_site.split(":", 1)[0]
        
    @property
    def site(self):
        return self.gene_site.split(":", 1)[1]

def main(count_files):
    d = dict()
    for cf in count_files:
        sample = os.path.basename(cf).split(".counts.txt.gz")[0]
        print >>sys.stderr, ">> processing", sample
        d[sample] = dict()
        for gene, cr_group in groupby(reader(cf, header=CountRecord), key=lambda x: x.gene):
            primary_site = None
            primary_site_count = 0
            # find the max among the group
            for cr in cr_group:
                if cr.count > primary_site_count:
                    primary_site = cr.site
                    primary_site_count = cr.count
            if primary_site:
                d[sample][primary_site] = True
    df = pandas.DataFrame(d)
    df.to_csv(sys.stdout, sep="\t", na_rep=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('count_files', nargs="+", help="count files per sample")
    args = p.parse_args()
    main(args.count_files)
