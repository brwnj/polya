#!/usr/bin/env python
# encoding: utf-8
"""
Create trackdb from list of folder contents.

+ &#43
- &#45
"""

import os
import re
import sys
import argparse
from string import Template
from toolshed import reader
from collections import defaultdict

COMPOSITE_TRACK_DEF = ("track PolyA\n"
                        "compositeTrack on\n"
                        "shortLabel Sites & Coverage\n"
                        "longLabel Sites and Coverage tracks\n"
                        "subGroup1 view Views PKS=Classified_Sites COV=Coverage SHT=Shifts\n"
                        "subGroup2 ptype PrimaryType UNK=Unknown CD34=CD34 CD71P=CD71&#43 CD71N=CD71&#45 K562=K562 TS=Tumor NBT=NormalTissue MCF7=MCF7 PK15=PK15 PK12=PK12 NA=NA\n"
                        "subGroup3 stype SecondaryType CD235P=CD235a&#43 CD235N=CD235a&#45 NOTX=No_Treatment HEMIN=Hemin PMA=PMA CONTROL=control ERP=ER&#43 ERN=ER&#45 PLACEBO=Placebo NA=NA\n"
                        "subGroup4 ttype TertiaryType PRP=PR&#43 PRN=PR&#45 TAM=Tamoxifen MPA=MPA NA=NA\n"
                        "subGroup5 qtype QuaternaryType HER2P=Her2&#43 HER2N=Her2&#45 NA=NA\n"
                        "subGroup6 strand Strand POS=Positive NEG=Negative U=Unstranded\n"
                        "subGroup7 inv Investigator MP=Pillai PK=Kabos\n"
                        "dimensions dimX=ptype dimY=stype dimA=inv dimB=strand dimC=ttype dimD=qtype\n"
                        "filterComposite dimA dimB dimC dimD\n"
                        "sortOrder view=-\n"
                        "type bed 6 .\n")

CLASSIFIED_PEAKS_VIEW = ("    track viewPeaks\n"
                            "    parent PolyA\n"
                            "    shortLabel Classified Sites\n"
                            "    view PKS\n"
                            "    visibility pack\n"
                            "    type bigBed 6 .\n")

COVERAGE_VIEW = ("    track viewCoverage\n"
                    "    parent PolyA\n"
                    "    shortLabel Coverage\n"
                    "    view COV\n"
                    "    visibility full\n"
                    "    maxHeightPixels 50:20:15\n"
                    "    type bigWig\n"
                    "    autoScale on\n"
                    "    alwaysZero on\n")

SITES_VIEW = """
track polyaSites
compositeTrack on
shortLabel Reference Sites
longLabel Classified Reference Sites
subGroup1 inv Investigator MP=Pillai PK=Kabos
visibility full
type bigBed 6 .

    track mp.test_sites
    parent polyaSites on
    subGroups inv=MP
    shortLabel Pillai 1,3,5
    longLabel Pillai: Class 1(a), 3(a), 5(a); exon model
    bigDataUrl MP.test_sites.bb
    color 215,48,39
    type bigBed 6 .

    track mp.all_sites
    parent polyaSites off
    subGroups inv=MP
    shortLabel Pillai All
    longLabel Pillai: Class 1(a), 2, 3(a), 4, 5(a), and 6; exon model
    bigDataUrl MP.all_sites.bb
    color 215,48,39
    type bigBed 6 .

    track mp.wholegene_sites
    parent polyaSites off
    subGroups inv=MP
    shortLabel WG:Pillai All
    longLabel Pillai: Class 1(a), 2, 3(a), 4, 5(a), and 6; whole gene model
    bigDataUrl MP.wholegene_sites.bb
    color 189,0,38
    type bigBed 6 .

    track pk.test_sites
    parent polyaSites on
    subGroups inv=PK
    shortLabel Kabos 1,3,5
    longLabel Kabos: Class 1(a), 3(a), 5(a); exon model
    bigDataUrl PK.test_sites.bb
    color 69,117,180
    type bigBed 6 .

    track pk.all_sites
    parent polyaSites off
    subGroups inv=PK
    shortLabel Kabos All
    longLabel Kabos: Class 1(a), 2, 3(a), 4, 5(a), and 6; exon model
    bigDataUrl PK.all_sites.bb
    color 69,117,180
    type bigBed 6 .

    track pk.wholegene_sites
    parent polyaSites off
    subGroups inv=PK
    shortLabel WG:Kabos All
    longLabel Kabos: Class 1(a), 2, 3(a), 4, 5(a), and 6; whole gene model
    bigDataUrl PK.wholegene_sites.bb
    color 37,52,148
    type bigBed 6 .
"""

# these will presumably have different sample subtypes
PK_FISHER_SHIFTS = ("track pkfishershifts\n"
                        "compositeTrack on\n"
                        "configurable on\n"
                        "shortLabel PK Fisher Shifts\n"
                        "longLabel Kabos: Observed Fisher shifts\n"
                        "subGroup1 strand Strand POS=Positive NEG=Negative\n"
                        "type bed 12 .\n")

MP_FISHER_SHIFTS = ("track mpfishershifts\n"
                        "compositeTrack on\n"
                        "configurable on\n"
                        "shortLabel MP Fisher Shifts\n"
                        "longLabel Pillai: Observed Fisher shifts\n"
                        "subGroup1 strand Strand POS=Positive NEG=Negative\n"
                        "type bed 12 .\n")

SITE_TEMPLATE = Template("        track $tname\n"
                            "        parent viewPeaks\n"
                            "        subGroups ptype=$ptype stype=$stype ttype=$ttype qtype=$qtype view=PKS inv=$inv strand=U\n"
                            "        shortLabel $slbl\n"
                            "        longLabel $llbl\n"
                            "        bigDataUrl $filename\n"
                            "        color $color\n"
                            "        type bigBed 6 .\n")

COVERAGE_TEMPLATE = Template("        track $tname\n"
                                "        parent viewCoverage\n"
                                "        subGroups ptype=$ptype stype=$stype ttype=$ttype qtype=$qtype view=COV inv=$inv strand=$strand\n"
                                "        shortLabel $slbl\n"
                                "        longLabel $llbl\n"
                                "        bigDataUrl $filename\n"
                                "        color $color\n"
                                "        type bigWig\n")

SHIFTS_TEMPLATE = Template("    track $tname\n"
                                "    parent $parent\n"
                                "    subGroups strand=$strand\n"
                                "    shortLabel $slbl\n"
                                "    longLabel $llbl\n"
                                "    bigDataUrl $filename\n"
                                "    color 37,52,148\n"
                                "    thickDrawItem on\n"
                                "    type bigBed 12 .\n")

def gsample(fname):
    return re.split('\.|_', fname)[0]

def gtname(fname):
    return os.path.splitext(fname)[0]

def ginv(sname):
    return "PK" if "PK" in sname else "MP"

def gstrand(fname):
    return "NEG" if "neg" in fname else "POS"

def gslbl(tname):
    # handling per sample sites
    if "classified" in tname:
        return "{sample}cpeaks".format(sample=tname.split("_")[0])
    # handling shifts
    if "_to_" in tname:
        tested = re.findall("(\d+)", tname)
        assert len(tested) == 2
        return "{first}>>{second}{strand}".format(first=tested[0], second=tested[1], strand="+" if "pos" in tname else "-")
    return tname[0:16]

def gllbl(tname, md):
    if "_to_" in tname:
        tested = re.findall("[MPK]+\d+", tname)
        assert len(tested) == 2
        return "{first} to {second} ({test}:{strand})".format(first=md[tested[0]]['translation'],
                                    second=md[tested[1]]['translation'],
                                    test="Fisher",
                                    strand=gstrand(tname))
    return "{translation} {strand}".format(translation=md[gsample(tname)]['translation'], strand=gstrand(tname))

def gshiftsparent(fname):
    if fname.startswith("PK"):
        return "pkfishershifts"
    if fname.startswith("MP"):
        return "mpfishershifts"

def flipstrand(fname):
    # since reads are being mapped to the opposite strand they belong to,
    # we are flipping the names, eg. 'pos' to 'neg'.
    tname = os.path.splitext(fname)[0]
    if "neg" in tname:
        tname = tname.replace("neg", "pos")
    elif "pos" in tname:
        tname = tname.replace("pos", "neg")
    return tname

def print_shifts(lst, pi, meta):
    """
    lst = list of file names
    pi = name of investigator or "pooled"
    meta = metadata as dict
    """
    for s in lst:
        tname = flipstrand(s)
        if pi == "pooled":
            if s.startswith("PK") or s.startswith("MP"): continue
            print SHIFTS_TEMPLATE.substitute(tname=tname,
                                    parent="pkfishershifts",
                                    strand=gstrand(tname),
                                    slbl=tname.split(".")[0],
                                    llbl=tname,
                                    filename=s)
        elif s.startswith(pi):
            print SHIFTS_TEMPLATE.substitute(tname=tname,
                                    parent=gshiftsparent(s),
                                    # stype="UNK",
                                    strand=gstrand(tname),
                                    slbl=gslbl(tname),
                                    llbl=gllbl(tname, meta),
                                    filename=s)

def main(folder, meta):
    filelist = os.listdir(folder)
    # group files by track type
    files = defaultdict(list)
    # metadata
    md = {}
    for l in reader(meta):
        md[l['alias']] = l

    for f in filelist:
        try:
            if md[gsample(f)]['exclude'] == "TRUE": continue
        except KeyError:
            # some files like trackDb will not be a valid key
            pass
        if f.endswith("bw"):
            files['coverage'].append(f)
            continue
        if "fisher" in f:
            files['fisher'].append(f)
            continue
        if "classified" in f:
            files['sites'].append(f)
            continue

    # composite track with classified peaks, coverage, and merged sites
    print COMPOSITE_TRACK_DEF

    # classified peaks for peter and manoj
    print CLASSIFIED_PEAKS_VIEW
    for s in files['sites']:
        sample = gsample(s)
        tname = gtname(s)
        print SITE_TEMPLATE.substitute(tname=tname,
                                    ptype=md[sample]['primary_type'],
                                    stype=md[sample]['secondary_type'],
                                    ttype=md[sample]['tertiary_type'],
                                    qtype=md[sample]['quaternary_type'],
                                    inv=ginv(s),
                                    slbl=gslbl(tname),
                                    llbl=gllbl(tname, md),
                                    filename=s,
                                    color=md[sample]['color'])
    print COVERAGE_VIEW
    # the bigwigs
    for s in files['coverage']:
        sample = gsample(s)
        tname = flipstrand(s)
        if s.startswith("PK") or s.startswith("MP"):
            print COVERAGE_TEMPLATE.substitute(tname=tname,
                                    ptype=md[sample]['primary_type'],
                                    stype=md[sample]['secondary_type'],
                                    ttype=md[sample]['tertiary_type'],
                                    qtype=md[sample]['quaternary_type'],
                                    inv=ginv(s),
                                    strand=gstrand(tname),
                                    slbl=gslbl(tname),
                                    llbl=gllbl(tname, md),
                                    filename=s,
                                    color=md[sample]['color'])
        else:
            print COVERAGE_TEMPLATE.substitute(tname=tname,
                                    ptype="NA",
                                    stype="NA",
                                    ttype="NA",
                                    qtype="NA",
                                    inv="PK",
                                    strand=gstrand(tname),
                                    slbl=tname,
                                    llbl=tname,
                                    filename=s,
                                    color="37,52,148")

    # should not change often, but could automate eventually
    print SITES_VIEW
    # shifts
    print PK_FISHER_SHIFTS
    print_shifts(files['fisher'], "PK", md)
    print_shifts(files['fisher'], "pooled", md)
    print MP_FISHER_SHIFTS
    print_shifts(files['fisher'], "MP", md)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("folder", help="files to convert to track definitions")
    p.add_argument("meta", help="google docs with sample info as tsv")
    args = p.parse_args()
    main(args.folder, args.meta)
