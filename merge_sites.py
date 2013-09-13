#!/usr/bin/env python
# encoding: utf-8
"""
Find the consensus of peaks among samples -- Assumes sample name is first in the
file name, delimited by either "." or "_" from the rest of the file name.
"""
import os
import sys
import tempfile
import os.path as op
import subprocess as sp
from toolshed import reader, nopen
from itertools import groupby, ifilterfalse, count, izip

class AnnotatedPeak(object):
    __slots__ = ['chrom','start','stop','name']
    def __init__(self, args, null="."):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        observed = []
        # classes from each sample reside after the name field
        for peak in args[4:]:
            if peak == null: continue
            observed.append(peak.rsplit(':', 1)[-1])
        # add most abundant classification onto name; name must be unique
        self.name += ":" + max(set(observed), key=observed.count)

    def __repr__(self):
        return "AnnotatedPeak({chr}:{name})".format(chr=self.chrom, name=self.name)

    def __str__(self):
        return "\t".join([getattr(self, s) for s in self.__slots__])

class PeakBed(object):
    __slots__ = ['chrom','start','stop','name','score','strand']
    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)

    def __repr__(self):
        return "PeakBed({chr}:{name})".format(chr=self.chrom, name=self.name)

    def __str__(self):
        return "\t".join([getattr(self, s) for s in self.__slots__])
    
    @property
    def pclass(self):
        return self.name.split(".")[1].lstrip("c")

def cleanup(files):
    """remove the files of a list."""
    for f in files:
        if type(f) is dict:
            for i in f.values():
                os.remove(i)
        else:
            os.remove(f)

def file_len(fname):
    p = sp.Popen(['wc', '-l', fname], stdout=sp.PIPE, stderr=sp.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def map_peak_class(annotated_peaks, merged_tmp):
    """deletes incoming temp file."""
    tmp = tempfile.mkstemp(suffix=".bed")[1]
    cmd = ("bedtools map -c 4 -o collapse -a {merged_tmp} "
            "-b {annotated_peaks} > {tmp}").format(**locals())
    sp.call(cmd, shell=True)
    # if no overlap exists, nothing is output by bedtools map
    if file_len(tmp) < 1:
        cleanup([tmp])
        return merged_tmp
    else:
        cleanup([merged_tmp])
        return tmp

def multi_intersect(files, cutoff):
    """files = {sample_name:file_path}"""
    sitestmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
    snames = [op.basename(f).split(".")[0].split("_")[0] for f in files]
    cmd = ("|bedtools multiinter -cluster -header "
                "-names {names} -i {files}").format(names=" ".join(snames),
                                                    files=" ".join(files))
    # apply cutoff, name peaks
    for i, l in enumerate(reader(cmd, header=True)):
        if int(l['num']) < cutoff: continue
        print >>sitestmp, "\t".join([l['chrom'], l['start'], l['end'],
                                        "peak_{i}".format(i=i)])
    sitestmp.close()
    # annotate the merged sites by intersecting with all of the files
    classtmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
    annotated_peaks = sitestmp.name
    # pull out peak classes from input files
    for f in files:
        annotated_peaks = map_peak_class(f, annotated_peaks)
    for peak in reader(annotated_peaks, header=AnnotatedPeak):
        if peak.name is None: continue
        print >>classtmp, "{chrom}\t{start}\t{stop}\t{name}\n".format(
                                chrom=peak.chrom, start=peak.start,
                                stop=peak.stop, name=peak.name)
    classtmp.close()
    return classtmp.name

def xref_to_dict(fname):
    d = {}
    for l in reader(fname, header=["from", "to"]):
        d[l['from']] = l['to']
    return d

def lparser(line, cols):
    return dict(zip(cols, line.strip().split("\t")))

def ret_item(line, cols, item):
    assert item in cols
    d = lparser(line, cols)
    return d[item]

def grouper(fp, cols):
    """yields group by gene name."""
    for k, g in groupby(fp, key=lambda t: ret_item(t, cols, "gene")):
        yield g

def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def get_out(l, n, xref):
    out = [l['chrom'],l['start'],l['stop']]
    gene = l['gene']
    # attempt to cross-reference an annotation
    # full, gene level translation messes with dexseq
    try:
        gene = "{gene}|{xref}".format(gene=l['gene'], xref=xref[l['gene']])
    except KeyError:
        # gene not in xref
        gene = l['gene']
    except TypeError:
        # xref not supplied
        gene = l['gene']
    out.append("p.c{pclass}.{gene}.{count}".format(
                                        pclass=l['peak'].rsplit(":", 1)[-1],
                                        gene=gene,
                                        count=n))
    out.extend(["0", l['strand']])
    return out

def add_slop(bed, sizes, n=0):
    """adds slop to bed entries. returns file name."""
    tmp = tempfile.mkstemp(suffix=".bed")[1]
    # slop output is not coordinate sorted
    cmd = ("bedtools slop -b {slop_amt} -i {bed}"
            " -g {sizes} | bedtools sort -i - > {out}").format(slop_amt=n,
                                                                bed=bed,
                                                                sizes=sizes,
                                                                out=tmp)
    sp.call(cmd, shell=True)
    return tmp

def intersect(ref, xref, peaks):
    if xref:
        xref = xref_to_dict(xref)
    # group the output by chr->gene->start
    cmd = ("|bedtools intersect -wb -a {peaks} -b {ref} "
                "| sort -k1,1 -k8,8 -k2,2n").format(**locals())
    cols = ['chrom','start','stop','peak','_chrom',
                '_start','_stop','gene','_score','strand']
    tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
    for g in grouper(nopen(cmd), cols):
        negs = []
        for i, l in enumerate(unique_everseen(\
                            g, lambda t: ret_item(t, cols, 'peak')), start=1):
            l = lparser(l, cols)
            # negative stranded sites
            if l['strand'] == "-":
                # need to count through them up, saving l each time
                negs.append(l)
                continue
            # positive stranded sites
            print >>tmp, "\t".join(get_out(l, i, xref))
        for i, l in izip(count(len(negs), -1), negs):
            print >>tmp, "\t".join(get_out(l, i, xref))
    tmp.close()
    return tmp.name

def filter_peaks(bed, classes):
    """removes peaks that do not match a class in classes."""
    classes = set(classes)
    for b in reader(bed, header=PeakBed):
        if b.pclass not in classes: continue
        print b

def main(ref, sizes, files, classes, xref, slop, cutoff):
    # get the overlapping peaks
    peak_regions = multi_intersect(files, cutoff)
    # add slop to gene model
    # ref_with_slop = add_slop(ref, sizes, slop)
    # annotate peaks with gene model
    # all_peaks = intersect(ref_with_slop, xref, peak_regions)
    all_peaks = intersect(ref, xref, peak_regions)
    # filter peaks by class
    filter_peaks(all_peaks, classes)
    # remove leftover temp files
    # cleanup([peak_regions, ref_with_slop, all_peaks])
    cleanup([peak_regions, all_peaks])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("ref", help="bed of regions with unique symbol (refseq, \
                                    ucsc, etc.) as name. it's important that \
                                    it's unique across transcripts if they're \
                                    present in the reference")
    # p.add_argument("sizes", help="chromosome sizes")
    p.add_argument("files", nargs="+", help="classified peaks")

    psites = p.add_argument_group("poly(A) sites")
    psites.add_argument("-c", metavar="CLASS", dest="classes", action="append",
            type=str, default=["1","1a"], choices=['1','1a','2','3','3a','4'],
            help="class of peaks used to generate consensus")
    psites.add_argument("-x", metavar="XREF", dest="xref", default=None,
            help="use xref value in place of ref name in polya site name; \
                    allowing rudimentary annotations. 2 columns: first is ref \
                    name, second is new value")

    pinter = p.add_argument_group("intersecting")
    # pinter.add_argument("-s", dest="slop", type=int, default=0,
    #         help="number of bases to add onto up and downstream ends of the \
    #                 gene model")
    pinter.add_argument("-n", dest="cutoff", type=int, default=2,
            help="number of samples containing called peak")

    args = p.parse_args()
    # revision: taking out slop option
    # main(args.ref, args.sizes, args.files, args.classes, args.xref, args.slop, args.cutoff)
    main(args.ref, args.files, args.classes, args.xref, args.cutoff)
