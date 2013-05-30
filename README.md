#Poly(A)

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
