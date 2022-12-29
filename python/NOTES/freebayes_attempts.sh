
## ------------------------ Parent bam files ---------------------------- ##


/analysis/CROTH/BAM/PMY2649_S112-rg.bam
/analysis/CROTH/BAM/PMY2650_S113-rg.bam


## ------------------------- Error messages from samtool index ------------------------------------- ##


/analysis/CROTH/BAM/PMY2557_S194-rg.bam
[E::hts_open_format] Failed to open file "/analysis/CROTH/BAM/PMY2557_S194-rg.bam" : Exec format error
samtools index: failed to open "/analysis/CROTH/BAM/PMY2557_S194-rg.bam": Exec format error
[E::hts_hopen] Failed to open file /analysis/CROTH/BAM/PMY2601_S299-rg.bam
[E::hts_open_format] Failed to open file "/analysis/CROTH/BAM/PMY2601_S299-rg.bam" : Exec format error
samtools index: failed to open "/analysis/CROTH/BAM/PMY2601_S299-rg.bam": Exec format error
[E::hts_hopen] Failed to open file /analysis/CROTH/BAM/PMY2701_S333-rg.bam
[E::hts_open_format] Failed to open file "/analysis/CROTH/BAM/PMY2701_S333-rg.bam" : Exec format error
samtools index: failed to open "/analysis/CROTH/BAM/PMY2701_S333-rg.bam": Exec format error
[E::hts_hopen] Failed to open file /analysis/CROTH/BAM/PMY2801_S45-rg.bam
[E::hts_open_format] Failed to open file "/analysis/CROTH/BAM/PMY2801_S45-rg.bam" : Exec format error
samtools index: failed to open "/analysis/CROTH/BAM/PMY2801_S45-rg.bam": Exec format error
[E::hts_hopen] Failed to open file /analysis/CROTH/BAM/PMY2901_S287-rg.bam
[E::hts_open_format] Failed to open file "/analysis/CROTH/BAM/PMY2901_S287-rg.bam" : Exec format error
samtools index: failed to open "/analysis/CROTH/BAM/PMY2901_S287-rg.bam": Exec format error


## ------------------------- Contig name of chromosome 13 ------------------------------ ##

CP003832.1 	

## -------------------------- Contig name of the mitochondria ------------------------- ## 

CP003834.1 

## ------------------- freebayes on genome for just parents ---------------------------------- ##

freebayes -p 1 -f /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta -L /analysis/CROTH/parentbams.txt | gzip > /analysis/CROTH/VCF/Bt22xFtc555.vcf.gz

## ------------------- freebayes on chromosome 1 for just parents ---------------------------------- ##

freebayes -p 1 -f /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta -t /analysis/CROTH/Bt22xFtc555-1/BED/CP003820.1.bed -Z -L /analysis/CROTH/parentbams.txt | gzip > /analysis/CROTH/VCF/Bt22xFtc555-1_CP003820.1.vcf.gz


## ------------------ Search linex system for a file -------------------- ##

grep --include=\*.{xml,py} -Rl ./ -e "Jason"


grep --include=\*.ipynb -Rl ./ -e "Cdx_QTL_SNPs.csv" 2> /dev/null

XL280x431_amoeba_assay_chromosome_7_QTL_ZOOM.pdf


grep --include=\*.ipynb -Rl ./ -e "XL280x431_amoeba_assay_chromosome_7_QTL_ZOOM.pdf" 2> /dev/null


grep --include=\*.ipynb -Rl ./ -e "KN99_knockouts.png" 2> /dev/null


grep --include=\*.ipynb -Rl ./ -e "H2O2_data.csv" 2> /dev/null


grep --include=\*.ipynb -Rl ./ -e "H2O2_scores.csv" 2> /dev/null


