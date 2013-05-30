#Poly(A)

Classifying peaks:
```
pclassify.py macs_peaks pos_bg neg_bg ref_fasta chrom_sizes \
    | bedtools sort -i - > classified_peaks
```

Building the poly(A) reference from classified peaks:
```
python merge_sites.py -n2 -c3 refseq.exon.bed.gz *classified_peaks.bed > passing_sites.bed
```

Getting the counts in order to run DEXSeq:
```
python read_counts.py bedgraph polya_sites.bed hg18.sizes > counts
```

Running DEXSeq:
```
python run_dexseq.py A.pos.counts B.pos.counts run_dexseq.R
```
