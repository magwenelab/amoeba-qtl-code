#!/bin/bash
for fn in tsauters/seqs;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i salmon_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 