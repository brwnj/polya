#!/usr/bin/env python
# encoding: utf-8
"""
Classify called peaks into categories as described Wang et al.:
http://rnajournal.cshlp.org/content/19/3/413.long

1   has A[A,T]TAAA; NOT A-rich downstream from this cleavage site
1a  has A[A,T]TAAA; NOT A-rich downstream from this cleavage site; A stretch immediately downstream of site
2   has A[A,T]TAAA; with A-rich sequence downstream
3   lacks A[A,T]TAAA; no A-rich region downstream
3a  lacks A[A,T]TAAA; no A-rich region downstream; A stretch immediately downstream of site
4   only downstream A-rich sequence

A gene on the '-' strand will have its peak called on '+' stranded reads. That
symbol should be accounted for when merging peaks from '+' and '-'. A bedgraph
for each '+' and '-' are properly assigned to corresponding peaks in this
script ('+' counts going to '-' peaks and vice versa).
"""
import os
import sys
import tempfile
import itertools
import subprocess as sp
from toolshed import reader

__version__ = "0.2"

def process_counts(pos, neg):
    """merges 2 bedgraphs, converts to bed with count in score field, and
    finally sorts. returns file name."""
    unsorted_tmp = tempfile.mkstemp(suffix=".txt")[1]
    tmp = tempfile.mkstemp(suffix=".bed")[1]
    # positive stranded counts belong to negative strand genes
    # add a negative symbol to mark for reverse complement later
    cmd = "gunzip -c %s" % pos if pos.endswith(".gz") else "cat %s" % pos
    cmd += """ | awk 'BEGIN{OFS=FS="\\t"}{print $1,$2,$3,".",$4,"-"}' > %s"""\
             % unsorted_tmp
    sp.call(cmd, shell=True)
    # negative stranded counts belong to positive strand genes
    cmd = "gunzip -c %s" % neg if neg.endswith(".gz") else "cat %s" % neg
    cmd += """ | awk 'BEGIN{OFS=FS="\\t"}{print $1,$2,$3,".",$4,"+"}' >> %s"""\
            % unsorted_tmp
    sp.call(cmd, shell=True)
    # sort the new bed file
    cmd = "bedtools sort -i %s > %s" % (unsorted_tmp, tmp)
    sp.call(cmd, shell=True)
    return tmp

def summit_start(starts, counts):
    """finds max position of starts given counts.
    
    >>> summit_start("20217423,20217424,20217425,20217426", "1,2,33,18")
    (20217425, 33)
    """
    starts = map(int, starts.split(","))
    counts = map(int, counts.split(","))
    summit_height = max(counts)
    return starts[counts.index(summit_height)], summit_height

def get_summits(bed, bedgraph):
    """overlap starts and counts with the peak regions and find the summit.
    """
    tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
    cmd = "|bedtools map -s -c 2 -o collapse -a %s -b %s |\
                 bedtools map -s -c 5 -o collapse -a - -b %s" % \
                 (bed, bedgraph, bedgraph)
    res = ["chr", "start", "stop", "name", "q", "strand", "starts", "counts"]
    for l in reader(cmd, header=res):
        # peak shifting isn't perfect so you'll get some peaks
        # where no reads overlap. those peaks typically have very few reads, so
        # just filter them out.
        if l['starts'] == ".":
            continue
        start, height = summit_start(l['starts'], l['counts'])
        fields = [l['chr'], start, start + 1, l['name'], height, l['strand']]
        tmp.write("\t".join(map(str, fields)) + "\n")
    tmp.close()
    return tmp.name

def add_slop(bed, sizes, slop):
    """adds slop to bed entries. returns file name."""
    tmp = tempfile.mkstemp(suffix=".bed")[1]
    cmd = "bedtools slop -b %d -i %s -g %s > %s" % (slop, bed, sizes, tmp)
    sp.call(cmd, shell=True)
    return tmp

def get_seqs(fasta, bed):
    """retrieves sequence in tab delimited file. returns file name."""
    tmp = tempfile.mkstemp(suffix=".txt")[1]
    # -s to reverse complement antisense
    cmd = "bedtools getfasta -s -tab -name -fi %s -bed %s -fo %s" % \
            (fasta, bed, tmp)
    sp.call(cmd, shell=True)
    return tmp

def seqs_to_dict(seqs):
    """name:sequence. returns dict."""
    d = {}
    for l in reader(seqs, header=["name","seq"]):
        d[l['name']] = l['seq']
    return d

def append_seqs(peaks, seqs):
    """takes tab delim seqs and appends seq column onto peak regions. output is
    written to a file. returns file name.
    """
    tmp = open(tempfile.mkstemp(suffix=".txt")[1], 'w')
    seqs = seqs_to_dict(seqs)
    res = ["chr", "start", "stop", "name", "count", "strand"]
    for l in reader(peaks, header=res):
        seq = seqs[l['name']].upper()
        fields = [l[i] for i in res] + [seq]
        tmp.write("\t".join(fields) + "\n")
    return tmp.name

def peak_category(seq, a_region_size, a_ratio, a_stretch, can_region):
    """determine the category of the peak. returns category as string.
    
    >>> peak_category("TTAAATAAATCTTCCTTTTAATTACTGGAAAAAATCTATTTT", 10, .65, 4, "-10,-30")
    '1'
    >>> peak_category("TTAAATAAATCTTCCTTCTTTTAAAAGGTTTACATCTATTTT", 10, .65, 4, "-10,-30")
    '1a'
    >>> peak_category("TTAAATAAATCTTCCTTTTAAAAAAAGGAAAAAATCTATTTT", 10, .65, 4, "-10,-30")
    '2'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTATCTGGAAGCAATCTATTTT", 10, .65, 4, "-10,-30")
    '3'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTAAAAGGAAGCAATCTATTTT", 10, .65, 4, "-10,-30")
    '3a'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAAAAAAAGGAAAAAATCTATTTT", 10, .65, 4, "-10,-30")
    '4'
    """
    canonical = ['AATAAA', 'ATTAAA']
    half_len = len(seq) / 2
    up = seq[:half_len]
    # these flip for string selection syntax
    range_right, range_left = map(int, can_region.split(","))
    # paper states PAS should be located in -10 to -30 region
    canonical_region = up[range_left:range_right]
    down = seq[half_len + 1:]
    found_pas = False
    found_arich = False
    alpha = False
    # no more than three canonical PAS located in the upstream window
    pas_count = 0
    for pas in canonical:
        pas_count += up.count(pas)
    if pas_count > 3:
            return False
    # has A[A,T]TAAA
    for pas in canonical:
        if canonical_region.find(pas) != -1:
            found_pas = True
            break
    # A-rich sequence downstream
    a_counts = down[:a_region_size].count("A")
    if a_counts > 0:
        a_content = a_counts / float(a_region_size)
        if a_content > a_ratio:
            found_arich = True
    # observed stretch of (A)s immediately downstream of site
    if down[:a_stretch].count("A") == a_stretch:
        alpha = True
    # return the category
    if found_pas and not found_arich:
        return "1a" if alpha else "1"
    if found_pas and found_arich:
        return "2"
    if not found_pas and not found_arich:
        return "3a" if alpha else "3"
    if not found_pas and found_arich:
        return "4"

def classify_peaks(bed, min_count, a_region_size, a_ratio, a_stretch,
                    can_region, out, size):
    """classify the peak sequences."""
    res = ["chrom", "start", "stop", "name", "count", "strand", "seq"]
    for i, l in enumerate(reader(bed, header=res), start=1):
        # filter out peaks with low read support
        if int(l['count']) < min_count: continue
        cat = peak_category(l['seq'], a_region_size, a_ratio, a_stretch, can_region)
        # no more than three canonical PAS located in the upstream window
        if not cat: continue
        # only output the summit
        l['start'] = int(l['start']) + size
        l['stop'] = l['start'] + 1
        l['name'] = "peak_%d:%s" % (i, cat)
        fields = [l[i] for i in out]
        print "\t".join(map(str, fields))

def main(args):
    if args.verbose:
        print >>sys.stderr, ">> munging count data"
    counts_tmp = process_counts(args.posbedgraph, args.negbedgraph)
    if args.verbose:
        print >>sys.stderr, ">> locating summit per peak"    
    # may want this file in the future
    summit_tmp = get_summits(args.peaks, counts_tmp)
    os.remove(counts_tmp)
    if args.verbose:
        print >>sys.stderr, ">> adding slop to summits"
    slop_tmp = add_slop(summit_tmp, args.sizes, args.wsize)
    os.remove(summit_tmp)
    if args.verbose:
        print >>sys.stderr, ">> retrieving peak sequences"
    # writes to a file; field 1 is 'chr1:1673768-1673869'
    seqs_tmp = get_seqs(args.fasta, slop_tmp)
    # appends the sequence column onto slop bed
    peak_tmp = append_seqs(slop_tmp, seqs_tmp)
    os.remove(slop_tmp)
    os.remove(seqs_tmp)
    if args.verbose:
        print >>sys.stderr, ">> classifying peaks"
    classify_peaks(peak_tmp, args.min_count, args.a_region, args.a_ratio,
                    args.a_stretch, args.can_region, args.out_cols, args.wsize)
    os.remove(peak_tmp)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('peaks', metavar="PEAKS", help="bed file of peak regions")
    p.add_argument('posbedgraph', metavar="POS_BEDGRAPH",
            help="bedgraph of counts from the positive strand")
    p.add_argument('negbedgraph', metavar="NEG_BEDGRAPH",
            help="bedgraph of counts from the negative strand")
    p.add_argument('fasta', metavar="FASTA", help="genomic reference sequence")
    p.add_argument('sizes', metavar="SIZES", help="genomic chromosome sizes")

    pslop = p.add_argument_group("slop settings")
    pslop.add_argument('-w', '--window-size', dest="wsize", type=int, default=50,
            help="number of nucleotides to add up- and \
                    downstream of peak site [%(default)s]")

    pclass = p.add_argument_group("peak classification")
    pclass.add_argument('--a-region-size', dest='a_region', type=int, default=10,
            help="number of bases to use when determining \
                    downstream A-content [%(default)s]")
    pclass.add_argument('--a-ratio-cutoff', dest='a_ratio', type=float,
            default=0.65, help="content ratio for A-region to be \
                                    considered A-rich [%(default)s]")
    pclass.add_argument('--a-stretch-length', dest='a_stretch', type=int,
            default=4, help="number of consecutive bases downstream of peak \
                                site to use when flagging classifications as \
                                alpha [%(default)s]")
    pclass.add_argument('--canonical-region', dest="can_region",
            default="-10,-30", help="narrowed upstream region in which \
                    canonical PAS should contained [%(default)s]")
    pclass.add_argument('--min-counts', dest="min_count", type=int, default=10,
            help="required read support for valid peak [%(default)s]")
    pclass.add_argument('-o', dest="out_cols", nargs="+",
            choices=["chrom", "start", "stop", "name", "count", "strand", "seq"],
            default=["chrom", "start", "stop", "name", "count", "strand"],
            help="define what is printed in the output; \
                    classification is incorporated into 'name' [BED6]")

    p.add_argument('--verbose', action='store_true',
            help="updates progress [%(default)s]")
    args = p.parse_args()
    main(args)
