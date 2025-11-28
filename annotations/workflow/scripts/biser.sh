#!/bin/bash
set -e

if [ -z $1 ]; then
        echo "Usage: ./biser.sh <BED> <FASTA> <out> "
        echo
        echo -e "\t<BED>:\tRepeatMasker output BED file to be used for softmasking."
        echo -e "\t<FASTA>:\tFASTA file with genome to annotate."
	echo -e "\t<threads>:\t number of threads to use for BISER."
        echo -e "\t<out>:\tprefix of the output files"
        echo
        echo "Softmask genome and run BISER to annotate SegDups."
        exit -1
fi


bed=$1
fasta=$2
threads=$3
outname=$4

echo "Softmask FASTA ..."
bedtools maskfasta -fi $fasta -bed $bed -fo ${outname}_softmasked.fa -soft &> ${outname}_softmasked.log
samtools faidx ${outname}_softmasked.fa

echo "Run BISER ..."
biser -o ${outname}_biser_segdups.bedpe -t ${threads} ${outname}_softmasked.fa --gc-heap 8G &> ${outname}_biser.log

echo "Convert to BED ..."
cat ${outname}_biser_segdups.bedpe | python3 workflow/scripts/biser-to-bed.py | bedtools sort | bedtools merge | python3 workflow/scripts/add-col.py > ${outname}_biser_segdups.bed

echo "Done."
