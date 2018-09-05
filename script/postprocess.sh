#! /bin/bash

# if [ ! -d ../result/matome ]
# then
#     mkdir -p ../result/matome 
# fi


python matome_omega.py ../output/kchiba_raw ../output/kchiba_merge/pcawg.savnet.result.txt ../db/CancerGeneSummary.proc.txt ../db/HGMD.CS.bed.gz 3 0.10



