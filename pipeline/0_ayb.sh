#!/usr/bin/env bash
#BSUB -J ayb_main
#BSUB -e ayb_main.%J.err
#BSUB -o ayb_main.%J.out
#BSUB -q normal
#BSUB -P pillai_kabos_polya

<<DOC
Barcodes file name format:
barcodes.LANE.txt, eg. barcodes.1.txt

Expected folder structure:
$HOME/projects/polya/data/20130305/Data/Intensities/L003/<cycles>/*.cif
DOC

set -o nounset -o errexit -o pipefail -x
source $HOME/projects/polya/bin/config.sh

lane=6
date=20140121

# 50 bp
bs=R6I14C37
# 100 bp
# bs=R6I14C86

cifs=$HOME/projects/polya/data/$date
cifs_dir=$cifs/Data/Intensities/L00${lane}/C1.1
fastq=$cifs/${date}_L00${lane}.fastq.gz
barcodes=$cifs/barcodes.${lane}.txt

jobids=""
for tile in `ls $cifs_dir | sed -rn 's/._._([0-9]+).cif/\1/p'`; do
    RUNSCRIPT=ayb.${tile}.$lane.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb_slave.${tile}.$lane" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e ayb.%J.$lane.err" >> $RUNSCRIPT
    echo "#BSUB -o ayb.%J.$lane.out" >> $RUNSCRIPT
    echo "#BSUB -n 8" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "#BSUB -P pillai_kabos_polya" >> $RUNSCRIPT
    echo "

AYB -b $bs -d cif -l debug -o $cifs -p 8 -i $cifs -r L${lane}T${tile} --format fastq
" >> $RUNSCRIPT

    job=$(bsub < $RUNSCRIPT)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
done

python -m bsub $jobids
cat $cifs/*${lane}*.fastq | gzip -c > $fastq
rm $cifs/?_${lane}_*.fastq
rm $cifs/ayb*.tab
fastq-multx -B $barcodes -m 2 -e $fastq -o $DATA/%.fq.gz
