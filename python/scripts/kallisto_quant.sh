#/bin/bash
#
#concatenate_sequence_multiple_lanes together
#
#
scriptargs="strain sequence_directory index.idx outdir #boots"
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
boots=$5

#make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

#find ${location} -name ${strain}*L001*R1*.fastq.gz -print

R1=$(find ${location} -name ${strain}_R1.fastq.gz)
R2=$(find ${location} -name ${strain}_R2.fastq.gz)


kallisto quant -i ${index} -o ${outdir}/${strain} -b ${boots} <(zcat ${R1}) <(zcat ${R2})

#echo "$LANE1_R1 $LANE1_R2 $LANE2_R1 $LANE2_R2"