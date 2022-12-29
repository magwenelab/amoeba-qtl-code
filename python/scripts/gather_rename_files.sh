#!/bin/bash
#
#concatenate_sequence_multiple_lanes.sh
#using subseq
#Requires indexed reference fasta and gff
#
#
scriptargs="strain home_directory outdir search_term(e.g *.vcf)"
nexpargs=4
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
outdir=$3
search=$4

#make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

for s in $(cat $strain) #need to do this in order to explode the samples supplied in a txt file
do
        find $location/$s -type f -name $search -print | xargs cp -t ../$outdir | xargs mv ../$outdir/$search ../$outdir/${s}_${search}
done