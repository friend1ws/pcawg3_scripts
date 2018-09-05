#! /usr/bin/env bash

for mutfile in `ls ../data/mutation_proc/*/*.mutation.txt`
do

    bfile=`basename $mutfile`
    sname=${bfile%.mutation.txt}

    dir1=`dirname $mutfile`
    ctype=`basename $dir1`
 
    echo -e "${ctype}\t${sname}"
        
    if [ ! -d ../output/intron_mut_annot/${ctype} ]; 
    then
        mkdir -p ../output/intron_mut_annot/${ctype}
    fi

    qsub subscript_intron_mut_annot/annot_mut.sh \
        ${mutfile} \
        ../output/intron_mut_annot/${ctype}/${sname}.mutation.annot.txt


done

