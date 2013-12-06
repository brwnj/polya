#Poly(A)

Scripts accompanying a pipeline to identify alternate polyadenylation sites.

#Workflow

1. Call bases using alternate base caller. https://github.com/timmassingham/AYB2
2. Move UMI from sequence onto the read name. https://github.com/brwnj/umitools

```
umitools trim --verbose $unprocessed_fastq $UMI | gzip -c > $fastq
```

3. Align the reads. We used Novoalign.
4. Remove duplicate UMIs across genome.

```
umitools rmdup --verbose $umibam $bam $UMI
```

5. Call peaks on stranded bams. https://github.com/taoliu/MACS

```
samtools view -hb -f 0x10 $bam > $negbam
samtools view -hb -F 0x10 $bam > $posbam
macs2 callpeak -t $negbam -n $negout --keep-dup auto --nomodel -s 25 --extsize 5
```

6. Sort out the peaks. Peaks called on the negative reads are from positive
stranded genes, so add sign according to file name (which is incorporated in 
the peak name).

```
zcat $negout.gz $posout.gz \
    | awk 'BEGIN{OFS=FS="\t"}{
            split($4, peak_basename, "/");
            $4 = peak_basename[length(peak_basename)];
            if($4~"neg"){$6="+"}
            else if($4~"pos"){$6="-"}
            print}' \
    | bedtools sort -i - \
    | gzip -c > $peak
```

7. Classify the peaks using bedgraphs from 5' read counts.

```
python classify_peaks.py $peak $posbg $negbg $FASTA $SIZES \
    | bedtools sort -i - \
    | gzip -c > $classified
````

8. Merge sites from across samples using only sites you intend to test.

```
python merge_sites.py -n2 -c3 -c3a -c5 -c5a -x $XREF $EXONS $peaks \
    | bedtools sort -i - \
    | gzip -c > $sites135
```

9. Get counts for each sample across the merged sites.
10. Run Fisher tests for desired contrasts.

```
python fisher_test.py a_counts b_counts \
    | gzip -c > a_to_b_fisher_result.txt
```
