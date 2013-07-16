#!/usr/bin/env python
# encoding: utf-8
"""
Find the consensus of peaks among samples -- Assumes sample name is first in the
file name, delimited by either "." or "_" from the rest of the file name.
"""
import os, sys, tempfile
import os.path as op
import subprocess as sp
from itertools import groupby, ifilterfalse, count, izip
from toolshed import reader, nopen

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

def filter_peaks(files, classes):
    """removes peaks that do not match a class in classes."""
    tmps = {}
    classes = set(classes)
    for f in files:
        sample = op.basename(f).split(".")[0].split("_")[0]
        tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
        res = ["chrom", "start", "stop", "name", "score", "strand"]
        i = 0
        for l in reader(f, header=res):
            c = int(l['name'].split(":")[1])
            if c not in classes: continue
            i += 1
            tmp.write("\t".join(l[i] for i in res) + "\n")
        tmp.close()
        # ensure at least one peak exists for this sample
        if i < 1:
            os.remove(tmp.name)
            continue
        tmps[sample] = tmp.name
    return tmps

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

def multi_intersect(slopfiles, origfiles, cutoff):
    """files = {sample_name:file_path}"""
    sitestmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
    names = " ".join(slopfiles.keys())
    paths = " ".join(slopfiles.values())
    cmd = ("|bedtools multiinter -cluster -header "
            "-names {names} -i {paths}").format(**locals())
    # apply cutoff, name peaks
    for i, l in enumerate(reader(cmd, header=True)):
        if int(l['num']) < cutoff: continue
        sitestmp.write("\t".join([l['chrom'], l['start'], l['end']]) + \
                            "\tpeak_{i}\n".format(i=i))
    sitestmp.close()
    # annotate the merged sites by intersecting with all of the files
    classtmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
    annotated_peaks = sitestmp.name
    # pull out peak classes from input files
    for sample_file in origfiles:
        annotated_peaks = map_peak_class(sample_file, annotated_peaks)
    for peak in reader(annotated_peaks, header=AnnotatedPeak):
        if peak.name is None: continue
        classtmp.write("{chrom}\t{start}\t{stop}\t{name}\n".format(
                    chrom=peak.chrom, start=peak.start,
                    stop=peak.stop, name=peak.name))
    classtmp.close()
    return classtmp.name

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

def get_out(l, n):
    out = [l['chrom'],l['start'],l['stop']]
    out.append("p.c{pclass}.{gene}.{count}".format(\
                                        pclass=l['peak'].rsplit(":", 1)[-1],\
                                        gene=l['gene'],\
                                        count=n))
    out.extend(["0", l['strand']])
    return out

def intersect(exons, peaks):
    # group the output by chr->gene->start
    cmd = ("|bedtools intersect -wb -a {peaks} -b {exons} "
                "| sort -k1,1 -k8,8 -k2,2n").format(**locals())
    cols = ['chrom','start','stop','peak','_chrom',
                '_start','_stop','gene','_score','strand']
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
            print "\t".join(get_out(l, i))
        for i, l in izip(count(len(negs), -1), negs):
            print "\t".join(get_out(l, i))

def main(exons, files, classes, cutoff):
    ## TODO:
    # filtering peaks should be done last in order to preserve site ids across
    # preferred class tracks
    tmps = filter_peaks(files, classes)
    peak_regions = multi_intersect(tmps, files, cutoff)
    intersect(exons, peak_regions)
    cleanup([tmps, peak_regions])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("exons", help="bed of unique exons with gene symbol as name")
    p.add_argument("files", nargs="+", help="classified peaks")

    psites = p.add_argument_group("poly(A) sites")
    psites.add_argument("-c", metavar="CLASS", dest="classes", action="append",
            type=int, default=[1], choices=[1,2,3,4],
            help="class of peaks used to generate consensus")

    pinter = p.add_argument_group("intersecting")
    pinter.add_argument("-n", dest="cutoff", type=int, default=2,
            help="number of samples containing called peak")

    args = vars(p.parse_args())
    main(**args)
