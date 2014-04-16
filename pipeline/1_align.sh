#!/usr/bin/env bash
#BSUB -J align[1-74]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 12
#BSUB -P pillai_kabos_polya

<<DOC
Trim the UMI from the FASTQ, align trimmed reads using Novoalign suppressing
all reads that align more than once, then remove UMI duplicates from the
alignment.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

unprocessed_fastq=$DATA/$sample.fq.gz
fastq=$DATA/$sample.umi.fq.gz

# trim the UMI
if [[ ! -f $fastq ]]; then
    umitools trim --verbose $unprocessed_fastq $UMI | gzip -c > $fastq
fi

results=$RESULTS/$sample
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

umibam=$RESULTS/$sample/$sample.UMIs_not_removed.bam
bam=$results/$sample.bam
stats=$results/$sample.alignment.txt
umibed=$results/$sample.umi_trimming.bed

# align the reads
if [[ ! -f $umibam ]]; then
    # just map the reads; no iterative trimming
    novoalign -d $NOVOIDX -f $fastq -o SAM -n 50 -r None -c 12 -k \
        2> $stats \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 8G \
        > $umibam
    samtools index $umibam
    bam2bw.py -5 -v -p pillai_kabos_polya $umibam
fi

# process the UMIs
if [[ ! -f $bam ]]; then
    umitools rmdup $umibam $bam > $umibed
    samtools index $bam
    bam2bw.py -5 -v -p pillai_kabos_polya $bam
fi
