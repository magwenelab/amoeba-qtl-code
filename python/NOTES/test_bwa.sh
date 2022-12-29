


## Wrote to /analysis/CROTH/SCRIPTS/test_bwa.sh
## Started at 1:55 pm
## finished vary quickly ... 
bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2585_S291_L002_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2585_S291_L002_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2585_S291_L002-sm.bam



## one of the larger files 
## Started at 2:02 pm
bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2691_S250_L001_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2691_S250_L001_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2691_S250_L001-sm.bam


