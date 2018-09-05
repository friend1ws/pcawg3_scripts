#! /bin/bash

if [ ! -d ../output/intron_mut_annot ]
then
    mkdir ../output/intron_mut_annot
fi


annot_utils gene --grc ../output/intron_mut_annot/gene.bed.gz

bedtools merge -i ../output/intron_mut_annot/gene.bed.gz > ../output/intron_mut_annot/gene.merged.bed


annot_utils exon --grc ../output/intron_mut_annot/exon.bed.gz

bedtools merge -i ../output/intron_mut_annot/exon.bed.gz > ../output/intron_mut_annot/exon.merged.bed


bedtools subtract -a ../output/intron_mut_annot/gene.merged.bed -b ../output/intron_mut_annot/exon.merged.bed > ../output/intron_mut_annot/intron.bed


# create cancer gene bed
python subscript_intron_mut_annot/create_cancer_gene_bed.py \
    ../output/intron_mut_annot/gene.bed.gz \
    ../db/cancer_gene.txt | \
    sort -k1,1 -k2,2n -k3,3n > \
    ../output/intron_mut_annot/gene.cg.bed

# create repeat masker bed file
python subscript_intron_mut_annot/proc_rmsk.py \
    ../db/rmsk.txt.gz | \
    sort -k1,1 -k2,2n -k3,3n > \
    ../output/intron_mut_annot/rmsk.bed

bgzip -c ../db/rmsk.bed > ../db/rmsk.bed.gz
tabix -p bed ../db/rmsk.bed.gz



