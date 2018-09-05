#! /usr/bin/env bash

for sfile in `ls ../data/ICGC-platinum-simple-files/*.indels.simple`
do

    bfile=`basename $sfile`
    ctype=${bfile%.indels.simple}

    echo $ctype

    if [ ! -d ../data/mutation_proc/${ctype} ]; 
    then
        mkdir -p ../data/mutation_proc/${ctype}
    fi

    python subscript_organize_mutation/organize_mutation.py \
        ../data/ICGC-platinum-simple-files/${ctype}.simple.gz \
        ../data/ICGC-platinum-simple-files/${ctype}.indels.simple \
        ../data/release_may2016.v1.4.tsv \
        ../data/mutation_proc/${ctype}

done


