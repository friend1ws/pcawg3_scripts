#! /bin/bash

if [ ! -d ../output/kchiba_raw ]
then
    mkdir -p ../output/kchiba_raw
fi

if [ ! -d ../output/kchiba_merge ]
then
    mkdir -p ../output/kchiba_merge
fi

if [ ! -d ../output/fdr ]
then
    mkdir -p ../output/fdr
fi

if [ ! -d ../output/motif ]
then
    mkdir -p ../output/motif
fi

if [ ! -d ../output/table ]
then
    mkdir -p ../output/table
fi


cp /home/omega3/icgc_pan_sanvet/*.savnet.result.txt ../output/kchiba_raw/

python execute_savnet_merge.py  ../output/kchiba_raw ../output/kchiba_merge/pcawg.savnet.result.txt ../db/CancerGeneSummary.proc.txt ../db/HGMD.CS.bed.gz 3 0.10

python check_cancer_type_fdr.py ../output/kchiba_raw/ 0.1 > ../output/fdr/cancer_fdr.txt  

python summarize_snv_info.py ../output/kchiba_merge/pcawg.savnet.result.txt ../output/motif/pcawg.savnet.motif_summary.txt /home/w3varann/database/GRCh37/GRCh37.fa ../db/hg19_spidex.bed.gz

Rscript seqlogo_summary.R


python proc_savnet.py ../output/kchiba_merge/pcawg.savnet.result.txt ../output/table/stable.raw.txt


