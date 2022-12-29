rm /analysis/CROTH/BAM/PMY2557_S194-rg.bam
rm /analysis/CROTH/BAM/PMY2601_S299-rg.bam
rm /analysis/CROTH/BAM/PMY2701_S333-rg.bam
rm /analysis/CROTH/BAM/PMY2801_S45-rg.bam
rm /analysis/CROTH/BAM/PMY2901_S287-rg.bam
gunzip /analysis/CROTH/BAM/PMY2557_S194.bam.gz
gunzip /analysis/CROTH/BAM/PMY2601_S299.bam.gz
gunzip /analysis/CROTH/BAM/PMY2701_S333.bam.gz
gunzip /analysis/CROTH/BAM/PMY2801_S45.bam.gz
gunzip /analysis/CROTH/BAM/PMY2901_S287.bam.gz
/home/croth/bin/./bamaddrg -b /analysis/CROTH/BAM/PMY2557_S194.bam -s PMY2557 -r S194.1 > /analysis/CROTH/BAM/PMY2557_S194-rg.bam
gzip /analysis/CROTH/BAM/PMY2557_S194.bam

/home/croth/bin/./bamaddrg -b /analysis/CROTH/BAM/PMY2601_S299.bam -s PMY2601 -r S299.45 > /analysis/CROTH/BAM/PMY2601_S299-rg.bam
gzip /analysis/CROTH/BAM/PMY2601_S299.bam

/home/croth/bin/./bamaddrg -b /analysis/CROTH/BAM/PMY2701_S333.bam -s PMY2701 -r S333.145 > /analysis/CROTH/BAM/PMY2701_S333-rg.bam
gzip /analysis/CROTH/BAM/PMY2701_S333.bam

/home/croth/bin/./bamaddrg -b /analysis/CROTH/BAM/PMY2801_S45.bam -s PMY2801 -r S45.245 > /analysis/CROTH/BAM/PMY2801_S45-rg.bam
gzip /analysis/CROTH/BAM/PMY2801_S45.bam

/home/croth/bin/./bamaddrg -b /analysis/CROTH/BAM/PMY2901_S287.bam -s PMY2901 -r S287.344 > /analysis/CROTH/BAM/PMY2901_S287-rg.bam
gzip /analysis/CROTH/BAM/PMY2901_S287.bam

samtools index /analysis/CROTH/BAM/PMY2557_S194-rg.bam
samtools index /analysis/CROTH/BAM/PMY2601_S299-rg.bam
samtools index /analysis/CROTH/BAM/PMY2701_S333-rg.bam
samtools index /analysis/CROTH/BAM/PMY2801_S45-rg.bam
samtools index /analysis/CROTH/BAM/PMY2901_S287-rg.bam
