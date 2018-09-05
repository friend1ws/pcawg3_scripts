#! /bin/bash

python subscript_intron_mut_annot/proc_savnet_result.py \
    ../output/kchiba_merge/pcawg.savnet.result.txt > \
    ../output/intron_mut_annot/savnet.creatin_mutation.txt

bash subscript_intron_mut_annot/annot_mut.sh \
    ../output/intron_mut_annot/savnet.creatin_mutation.txt \
    ../output/intron_mut_annot/savnet.creatin_mutation.annot.txt

python subscript_intron_mut_annot/count_savnet_annot.py \
    ../output/intron_mut_annot/savnet.creatin_mutation.annot.txt > \
    ../output/intron_mut_annot/savnet.creatin_mutation.annot_summary.txt

Rscript subscript_intron_mut_annot/alt_junc_pos_log.R

Rscript subscript_intron_mut_annot/deep_intron_annot_plot.R


