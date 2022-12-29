bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2649_S112_L001_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2649_S112_L001_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2649_S112_L001-sm.bam

bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2649_S112_L002_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2649_S112_L002_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2649_S112_L002-sm.bam

bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2650_S113_L001_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2650_S113_L001_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2650_S113_L001-sm.bam

bwa mem -a -M /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2650_S113_L002_R1_001.fastq.gz /data/sequence-data/Crypto_Bt22xFtc555_Oct2020/PMY2650_S113_L002_R2_001.fastq.gz | samtools view -F 4 -b | samtools sort -o /analysis/CROTH/SAM/PMY2650_S113_L002-sm.bam

