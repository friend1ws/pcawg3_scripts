#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -e log -o log

INPUT=$1
OUTPUT=$2

cut -f 1-5 ${INPUT} | bedtools intersect -a - -b ../output/intron_mut_annot/gene.bed.gz -wb > ${OUTPUT}.tmp.gene.bed

cut -f 1-5 ${INPUT} | bedtools intersect -a - -b ../output/intron_mut_annot/intron.bed -wb > ${OUTPUT}.tmp.intron.bed

cut -f 1-5 ${INPUT} | bedtools intersect -a - -b ../output/intron_mut_annot/gene.cg.bed -wb > ${OUTPUT}.tmp.cg.bed

cut -f 1-5 ${INPUT} | bedtools intersect -a - -b ../output/intron_mut_annot/rmsk.bed -wb > ${OUTPUT}.tmp.rmsk.bed

python subscript_intron_mut_annot/add_annot.py \
    ${INPUT} ${OUTPUT}.tmp.gene.bed \
    ${OUTPUT}.tmp.intron.bed ${OUTPUT}.tmp.cg.bed \
    ${OUTPUT}.tmp.rmsk.bed > \
    ${OUTPUT} 


rm -rf ${OUTPUT}.tmp.gene.bed
rm -rf ${OUTPUT}.tmp.intron.bed
rm -rf ${OUTPUT}.tmp.cg.bed
rm -rf ${OUTPUT}.tmp.rmsk.bed



