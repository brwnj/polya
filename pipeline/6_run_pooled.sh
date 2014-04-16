#!/usr/bin/env bash
#BSUB -J run_pooled
#BSUB -e run_pooled.%J.err
#BSUB -o run_pooled.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
testing pooled.

files are output to $RESULT/pooled_results
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

cd $POOLEDRESULTS

# this script will filter UMI_not_removed bedgraphs
bedgraphs=$HOME/projects/polya/results/common/*/*.bedgraph.gz
python $HOME/projects/polya/bin/get_pooled_coverage.py -v $METADATA $bedgraphs
for bedgraph in *.bedgraph; do
    bedGraphToBigWig $bedgraph $SIZES ${bedgraph/.bedgraph/.bw}
    gzip -f $bedgraph
done

# get the counts for each pool from the bedgraphs
for bg in *.bedgraph.gz; do
    # pull strand from bedgraph file name
    strand="pos"
    if [[ "${bg#*pos}" == "$bg" ]]; then
        strand="neg"
    fi
    countsout=${bg/.bedgraph.gz}.counts.txt.gz
    slopsites=$POLYASITES/PK.test_sites.slop.$SLOP.$strand.bed.gz

    python $BIN/read_counts.py $bg $slopsites | gzip -c > $countsout
done

ext=counts.txt.gz
# run fisher tests across desired comparisons
for strand in pos neg; do
    bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.$strand.$ext TS-ERP.$strand.$ext | gzip -c > NBT_to_TS-ERP.$strand.fisher.txt.gz" &
    bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.$strand.$ext TS.$strand.$ext | gzip -c > NBT_to_TS.$strand.fisher.txt.gz" &
    bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.$strand.$ext TS-ERN.$strand.$ext | gzip -c > NBT_to_TS-ERN.$strand.fisher.txt.gz" &
    bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py TS-ERP.$strand.$ext TS-ERN.$strand.$ext | gzip -c > TS-ERP_to_TS-ERN.$strand.fisher.txt.gz" &
done
wait

# convert fisher test results to bed format
pksites=$POLYASITES/PK.test_sites.bed.gz
for strand in pos neg; do
    bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS.$strand.fisher.txt.gz $pksites | sort -k1,1 -k2,2n -k3,3n | uniq > NBT_to_TS.$strand.fisher.bed" &
    bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERP.$strand.fisher.txt.gz $pksites | sort -k1,1 -k2,2n -k3,3n | uniq > NBT_to_TS-ERP.$strand.fisher.bed" &
    bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERN.$strand.fisher.txt.gz $pksites | sort -k1,1 -k2,2n -k3,3n | uniq > NBT_to_TS-ERN.$strand.fisher.bed" &
    bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py TS-ERP_to_TS-ERN.$strand.fisher.txt.gz $pksites | sort -k1,1 -k2,2n -k3,3n | uniq > TS-ERP_to_TS-ERN.$strand.fisher.bed" &
done
wait

# convert bed to bb
for strand in pos neg; do
    bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $SIZES NBT_to_TS.$strand.fisher.bed" &
    bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $SIZES NBT_to_TS-ERP.$strand.fisher.bed" &
    bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $SIZES NBT_to_TS-ERN.$strand.fisher.bed" &
    bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $SIZES TS-ERP_to_TS-ERN.$strand.fisher.bed" &
done
wait

# copy bws and bbs over to the hub
cp *.bb $HUB/hg19
cp *.bw $HUB/hg19
