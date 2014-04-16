#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
pull all necessary files and convert as necessary
DOC

set -o nounset -x
source $HOME/projects/polya/bin/config.sh

cd $HUB/hg19
# copy over BWs
cp $RESULTS/*/*.bw .
rm *UMIs_not_removed*
# copy over classified peaks
cp $RESULTS/*/*classified.bed.gz .
gunzip -f *.classified.bed.gz
# convert to BB6
for bed in *.classified.bed; do bed2bb.py --type bed6 $SIZES $bed; done
# delete classified peaks.bed
rm *.classified.bed

# copy over sites beds
cp $POLYASITES/*sites*.bed.gz .
# convert to BB
gunzip -f *sites*.bed.gz
for bed in *sites*.bed; do bed2bb.py --type bed6 $SIZES $bed; done
# these file names are hard coded into the python script: generate_trackdb.py
rm *sites*.bed

for f in $FISHERRESULTS/*fisher*txt.gz; do
    fbase=$(basename $f)
    fisherbed=${f/.txt.gz/.bed}
    fisherbb=${f/.txt.gz/.bb}
    if [[ ! -f $fisherbb ]]; then
        # create beds of fisher test results
        python $BIN/visualize_fisher_shifts.py $f $POLYASITES/${fbase:0:2}.test_sites.bed.gz \
            | bedtools sort -i - > $fisherbed
        # convert bed to bigbed
        bed2bb.py --type bed12 $SIZES $fisherbed
        # delete the bed
        rm $fisherbed
    fi
    # cp over the bb
    cp $fisherbb $HUB/hg19/
done

# run generate_trackdb.py > trackDb.txt
python $HOME/projects/polya/bin/generate_trackdb.py . $METADATA > trackDb.txt
