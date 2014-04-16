#! /usr/bin/env bash
#BSUB -J counts[1-74]
#BSUB -e counts.%J.%I.err
#BSUB -o counts.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

<<DOC
get the counts. overwrite existing count files.
DOC

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULTS/$sample

for strand in pos neg; do
    bedg=$results/$sample.$strand.bedgraph.gz
    countsout=$results/$sample.$strand.counts.txt.gz
    slopsites=$POLYASITES/${sample:0:2}.test_sites.slop.$SLOP.$strand.bed.gz

    # write counts to a temp file
    python $BIN/read_counts.py $bedg $slopsites | gzip -c > $countsout.tmp.gz
    # filter duplicate gene entries from the counts
    # input is sorted by gene
    python $HOME/projects/polya/bin/filter_duplicate_counts.py $countsout.tmp.gz | gzip -c > $countsout
    # remove temp
    rm $countsout.tmp.gz
done

# python ~/devel/polya/misc/primary_site_annotation_table.py ~/projects/polya/results/common/*/*counts.txt.gz | gzip -c > primary_site_annotation_table.txt.gz
