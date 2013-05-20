#Poly(A)

Building the poly(A) reference from classified peaks:
```
python merge_sites.py -b 0 ref/refseq.exon.bed.gz ref/hg18.sizes *classified.bed.gz > sites.bed
```

Getting the counts in order to run DEXSeq:
```
bedtools map -s -c 4 -o max -null 0 -a sites_with_slop.bed -b MP55.pos.bedgraph.gz | cut -f 4,7
```

Generating "replicates":
```
for f in *.counts;\
    do awk 'BEGIN{OFS=FS="\t"}{foo=int(rand()*10);if($2!=0){print $1,$2+foo}else{print}}' $f\
    > ${f%.counts}.x.counts;\
done
```
