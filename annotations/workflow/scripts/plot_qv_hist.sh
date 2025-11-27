#!/bin/bash
set -e

if [ -z $1 ]; then
        echo "Usage: ./plot_qv_hist.sh <BED> <TSV> <out> "
        echo
        echo -e "\t<BED>:\tRepeatMasker output BED file."
        echo -e "\t<TSV>:\tWindow-based QVs"
        echo -e "\t<out>:\tprefix of the output files"
        echo
        echo "Plot QV histogram with repeat annotations."
        exit -1
fi


rb=$1
w=$2
outname=$3


# generate separate BED file for each repeatclass
cat $rb | python3 workflow/scripts/split-by-repeatclass.py $outname

# convert QV tsv to a BED file
cat $w | python3 workflow/scripts/convert-to-bed.py > ${outname}_windows.bed

# annotate QV BED
bedtools annotate -i ${outname}_windows.bed -files ${outname}_DNA_merged.bed ${outname}_LINE_merged.bed ${outname}_SINE_merged.bed ${outname}_LTR_merged.bed ${outname}_Satellite_merged.bed ${outname}_OTHER_merged.bed ${outname}_SegDup.bed -names DNA LINE SINE LTR Satellite OtherRepeats SegDups | bedtools sort -header > ${outname}_windows_annotated.bed

# plot
cat ${outname}_windows_annotated.bed | python3 workflow/scripts/plot_annotated_histogram.py ${outname}_histogram.pdf
