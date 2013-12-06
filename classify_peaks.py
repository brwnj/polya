#!/usr/bin/env python
# encoding: utf-8
"""
Classify called peaks into categories as described Wang et al.:
http://rnajournal.cshlp.org/content/19/3/413.long

1   has A[A,T]TAAA; not A-rich downstream
1a  has A[A,T]TAAA; not A-rich downstream; A-stretch immediately downstream
2   has A[A,T]TAAA; A-rich downstream sequence
3   lacks A[A,T]TAAA, non-canonical site; not A-rich downstream
3a  lacks A[A,T]TAAA, non-canonical site; not A-rich downstream; A-stretch immediately downstream
4   only downstream A-rich sequence

Additional sub-category to further characterize class 3 sites:

5   has non-canonical site; not A-rich downstream
5a  has non-canonical site; not A-rich downstream; A-stretch immediately downstream

A gene on the '-' strand will have its peak called on '+' stranded reads. That
symbol should be accounted for when merging peaks from '+' and '-'. A bedgraph
for each '+' and '-' are properly assigned to corresponding peaks in this
script ('+' counts going to '-' peaks and vice versa).
"""

import os
import re
import sys
import doctest
import tempfile
import itertools
import subprocess as sp
from toolshed import reader
from argparse import ArgumentParser, RawDescriptionHelpFormatter

os.environ['TERM'] = 'linux'
__version__ = "0.3"

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

def divide_peak_sequence(seq, canonical_region):
    """
    Split the peak sequence relative to its center, then slice out the defined
    canonical target region from the upstream sequence.

    >>> divide_peak_sequence("TTAAATAAATCTTCCTTTTAATTTACTGGAAAAAATCTATTT", [-5,-19])
    ('TTAAATAAATCTTCCTTTTAA', 'TTACTGGAAAAAATCTATTT', 'AAATAAATCTTCCT')
    """
    half_len = len(seq) / 2
    upstream = seq[:half_len]
    # flip canonical_region for string selection syntax
    range_right, range_left = canonical_region
    # paper states PAS should be located in -10 to -30 region
    # though we've relaxed that to -5 and -40 as our defaults
    target_region = upstream[range_left:range_right]
    downstream = seq[half_len + 1:]
    return upstream, downstream, target_region

def test_a_rich(seq, a_region_size, a_ratio):
    """
    Testing the number of As found in just downstream of the PAS.

    seq                 target sequence found downstream of the PAS
    a_region_size       number of bases to look at from the PAS
    a_ratio             A content >= to this ratio will be considered an A-stretch

    >>> test_a_rich("AAAAAGGGGGTTTTTTTTTT", 10, .5)
    True
    >>> test_a_rich("AAAAAGGGGGTTTTTTTTTT", 10, .6)
    False
    """
    found_a_rich = False
    # A-rich sequence downstream
    a_counts = seq[:a_region_size].count("A")
    if a_counts > 0:
        a_content = a_counts / float(a_region_size)
        if a_content >= a_ratio:
            found_a_rich = True
    return found_a_rich

def test_for_a_stretch(seq, a_stretch):
    """
    Test whether or not the immediate downstream bases of the PAS are A
    throughout the a_stretch distance.

    >>> test_for_a_stretch("AAAAATTTTTTTTTT", 5)
    'a'
    >>> test_for_a_stretch("AAAATTTTTTTTTTT", 5)
    ''
    """
    return 'a' if seq[:a_stretch].count("A") == a_stretch else ''

def total_canonical_hits(canonical_seqs, upstream_sequence):
    """
    Count the number of times canonical sequences appear in the upstream
    sequence.

    >>> total_canonical_hits(['AATAAA','ATTAAA'], "TTAAATAAATCTTTAAATCTTTAAAATCT")
    1
    >>> total_canonical_hits(['AATAAA','ATTAAA'], "TTAAATAAATTTAAATTAAATTAAATAAATAATAAACT")
    4
    """
    hits = 0
    for seq in canonical_seqs:
        hits += upstream_sequence.count(seq)
    return hits

def peak_category(seq, a_region_size, a_ratio, a_stretch, can_region, noncan_seqs):
    """
    Determine classification for each peak region. It first attempts to
    classify to a canonical PAS. Any resultant class 3 peak is then checked
    against the non-canonical sequence list.

    Class 2:
    No more than three canonical PAS should exist in the 50-nt window upstream
    of the cleavage site because multiple occurrences of PAS is likely due to
    interspersed non-A nucleotide in an A-rich sequence segment... indicative
    of internal priming events.

    >>> peak_category("TTAAATAAATCTTCCTTTTAATTACTGGAAAAAATCTATTTT", 10, .65, 4, [-10,-30], False)
    '1'
    >>> peak_category("TTAAATAAATCTTCCTTCTTTTAAAAGGTTTACATCTATTTT", 10, .65, 4, [-10,-30], False)
    '1a'
    >>> peak_category("TTAAATAAATCTTCCTTTTAAAAAAAGGAAAAAATCTATTTT", 10, .65, 4, [-10,-30], False)
    '2'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTATCTGGAAGCAATCTATTTT", 10, .65, 4, [-10,-30], False)
    '3'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTAAAAGGAAGCAATCTATTTT", 10, .65, 4, [-10,-30], False)
    '3a'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAAAAAAAGGAAAAAATCTATTTT", 10, .65, 4, [-10,-30], False)
    '4'

    class 2 sequence with more than 3 PAS
    >>> peak_category("AATAAAAATAAAAATAAAATTAAAAAAAAAAAAAAAAAAAAAAAAATT", 10, .65, 4, [-10,-30], False)
    ''

    >>> noncan_seqs = ['AAGAAA']
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTATCTGGAAGCAATCTATTTT", 10, .65, 4, [-10,-30], noncan_seqs)
    '5'
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTAAAAGGAAGCAATCTATTTT", 10, .65, 4, [-10,-30], noncan_seqs)
    '5a'

    sequence has non-canonical hit, but is A-rich downstream
    >>> peak_category("TTAAAGAAATCTTCCTTTTAGTAAAAAAAAAAAATCTATTTT", 10, .65, 4, [-10,-30], noncan_seqs)
    '4'
    >>>
    """
    site_class = ""

    # the canonical sequences used to determine class 1
    canonical = ['AATAAA', 'ATTAAA']

    # split up the sequence into logical portions
    upstream, downstream, target_region = divide_peak_sequence(seq, can_region)

    # analyze downstream sequence
    downstream_a_rich = test_a_rich(downstream, a_region_size, a_ratio)
    alpha_subcategory = '' if downstream_a_rich else test_for_a_stretch(downstream, a_stretch)

    # test the canonical sequences first
    for query in canonical:
        if target_region.find(query) != -1:
            site_class = "2" if downstream_a_rich else "1"
            break
        elif not downstream_a_rich:
            site_class = "3"
        else:
            site_class = "4"

    # need to test class 2s further: see docstring
    if site_class == "2":
        hits = total_canonical_hits(canonical, upstream)
        if hits > 3:
            site_class = ""

    # test non-canonical sequences against class 3 sites
    if site_class == "3" and noncan_seqs:
        for query in noncan_seqs:
            if target_region.find(query) != -1:
                site_class = "5"
                break

    return site_class + alpha_subcategory

def classify_peaks(bed, min_count, a_region_size, a_ratio, a_stretch,
                    can_region, noncan_seqs, out_cols, size):
    """classify the peak sequences."""
    res = ["chrom", "start", "stop", "name", "count", "strand", "seq"]
    for i, l in enumerate(reader(bed, header=res), start=1):
        # filter out peaks with low read support
        if int(l['count']) < min_count: continue
        cat = peak_category(l['seq'], a_region_size, a_ratio, a_stretch, can_region, noncan_seqs)
        # no more than three canonical PAS located in the upstream window
        if not cat: continue
        # only output the summit
        l['start'] = int(l['start']) + size
        l['stop'] = l['start'] + 1
        l['name'] = "peak_%d:%s" % (i, cat)
        fields = [l[i] for i in out_cols]
        print "\t".join(map(str, fields))

def main(args):
    print >>sys.stderr, ">> munging count data"
    counts_tmp = process_counts(args.posbedgraph, args.negbedgraph)

    print >>sys.stderr, ">> locating summit per peak"
    # may want this file in the future
    summit_tmp = get_summits(args.peaks, counts_tmp)
    os.remove(counts_tmp)

    print >>sys.stderr, ">> adding slop to summits"
    slop_tmp = add_slop(summit_tmp, args.sizes, args.wsize)
    os.remove(summit_tmp)

    print >>sys.stderr, ">> retrieving peak sequences"
    # writes to a file; field 1 is 'chr1:1673768-1673869'
    seqs_tmp = get_seqs(args.fasta, slop_tmp)
    # appends the sequence column onto slop bed
    peak_tmp = append_seqs(slop_tmp, seqs_tmp)
    os.remove(slop_tmp)
    os.remove(seqs_tmp)

    print >>sys.stderr, ">> classifying peaks"
    classify_peaks(peak_tmp, args.min_count, args.a_region, args.a_ratio, args.a_stretch, args.can_region, args.noncan_seqs, args.out_cols, args.wsize)
    os.remove(peak_tmp)

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, version=__version__, formatter_class=RawDescriptionHelpFormatter)

    p.add_argument('peaks', metavar="PEAKS", help="bed file of peak regions")
    p.add_argument('posbedgraph', metavar="POS_BEDGRAPH", help="bedgraph of counts from the positive strand")
    p.add_argument('negbedgraph', metavar="NEG_BEDGRAPH", help="bedgraph of counts from the negative strand")
    p.add_argument('fasta', metavar="FASTA", help="genomic reference sequence")
    p.add_argument('sizes', metavar="SIZES", help="genomic chromosome sizes")

    pslop = p.add_argument_group("slop settings")
    pslop.add_argument('-w', '--window-size', dest="wsize", type=int, default=50, help="number of nucleotides to add up- and downstream of peak site [%(default)s]")

    pclass = p.add_argument_group("peak classification")
    pclass.add_argument('--a-region-size', dest='a_region', type=int, default=10, help="number of bases to use when determining downstream A-content [%(default)s]")
    pclass.add_argument('--a-ratio-cutoff', dest='a_ratio', type=float, default=0.65, help="content ratio for A-region to be considered A-rich [%(default)s]")
    pclass.add_argument('--a-stretch-length', dest='a_stretch', type=int, default=4, help="number of consecutive bases downstream of peak site to use when flagging classifications as alpha [%(default)s]")
    pclass.add_argument('--canonical-region', dest="can_region", type=int, nargs=2, default=[-5,-40], help="narrowed upstream region in which canonical PAS should contained [%(default)s]")
    pclass.add_argument('--min-counts', dest="min_count", type=int, default=10, help="required read support for valid peak [%(default)s]")
    pclass.add_argument('--output-columns', dest="out_cols", nargs="+", choices=["chrom", "start", "stop", "name", "count", "strand", "seq"], default=["chrom", "start", "stop", "name", "count", "strand"], help="define what is printed in the output; classification is incorporated into 'name' [BED6]")
    pclass.add_argument('-s', '--noncanonical-seqs', dest="noncan_seqs", help="if specified, class 3 sites will be further subdivided using these 6 nt sequences, eg. -s AGTAAA,TATAAA,TTTAAA, giving class 5 peaks")

    args = p.parse_args()

    if not args.can_region[0] > args.can_region[1]:
        print >>sys.stderr, "Canonical Region:\nLarger number first, e.g. -5 -40"
        sys.exit(1)

    if args.noncan_seqs:
        r = re.compile(r"^[ACTGactg]{6}$")
        args.noncan_seqs = [seq.upper() for seq in args.noncan_seqs.split(",") if r.match(seq)]

    if doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE).failed == 0:
        main(args)
