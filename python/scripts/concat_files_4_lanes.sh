#/bin/bash
#
#concatenate_sequence_multiple_lanes together
#
#
scriptargs="strain sequence_directory outdir"
nexpargs=3
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

#make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

#find ${location} -name ${strain}*L001*R1*.fastq.gz -print

LANE1_R1=$(find ${location} -name ${strain}_*L001*R1*.fastq.gz)
LANE1_R2=$(find ${location} -name ${strain}_*L001*R2*.fastq.gz)
LANE2_R1=$(find ${location} -name ${strain}_*L002*R1*.fastq.gz)
LANE2_R2=$(find ${location} -name ${strain}_*L002*R2*.fastq.gz)
LANE3_R1=$(find ${location} -name ${strain}_*L003*R1*.fastq.gz)
LANE3_R2=$(find ${location} -name ${strain}_*L003*R2*.fastq.gz)
LANE4_R1=$(find ${location} -name ${strain}_*L004*R1*.fastq.gz)
LANE4_R2=$(find ${location} -name ${strain}_*L004*R2*.fastq.gz)

cat $LANE1_R1 $LANE2_R1 $LANE3_R1 $LANE4_R1 > ${outdir}/${strain}_R1.fastq.gz
cat $LANE1_R2 $LANE2_R2 $LANE3_R1 $LANE4_R2> ${outdir}/${strain}_R2.fastq.gz

#echo "$LANE1_R1 $LANE1_R2 $LANE2_R1 $LANE2_R2"