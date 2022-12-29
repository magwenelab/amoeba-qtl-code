#!/bin/bash

#/bin/bash
#
#my_salmon_quant
#
#
scriptargs="strain sequence_directory index_file outdir #cores"
nexpargs=5
nargs=$#

if [[ $nargs -ne $nexpargs ]]
then
    echo
    echo "Usage: $(basename $0) scriptargs"
    echo
    exit
fi

strain=$1
location=$2
index=$3
outdir=$4
cores=$5

#make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

salmon quant -i ${index} -l A \
         -1 ${location}/${strain}_R1.fastq.gz \
         -2 ${location}/${strain}_R2.fastq.gz \
         -p ${cores} --validateMappings -o ${outdir}/${strain}_quant