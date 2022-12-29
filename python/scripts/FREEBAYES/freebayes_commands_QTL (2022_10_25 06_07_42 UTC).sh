## ---------- Call variants on chromosome 8 (CP003827.1) QTL -------------- ##

freebayes -p 1 -f /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta -r CP003827.1:660000-770000 -Z -L /analysis/CROTH/Bt22xFtc555-1/SCRIPTS/qtllistofbams.txt | gzip > /analysis/CROTH/VCF/Bt22xFtc555-1_CP003827.1_Chr8_QTL.vcf.gz


## ---------- Call variants on chromosome 5 (CP003824.1) QTL -------------- ##

freebayes -p 1 -f /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta -r CP003824.1:1015000-1108000 -Z -L /analysis/CROTH/Bt22xFtc555-1/SCRIPTS/qtllistofbams.txt | gzip > /analysis/CROTH/VCF/Bt22xFtc555-1_CP003824.1_Chr5_QTL.vcf.gz

## ------------ Call variants just for 1 kb of Chromosome 8

freebayes -p 1 -f /analysis/CROTH/REF/FungiDB-46_CneoformansH99_Genome.fasta -r CP003827.1:660000-661000 -Z -L /analysis/CROTH/Bt22xFtc555-1/SCRIPTS/qtllistofbams.txt | gzip > /analysis/CROTH/VCF/Bt22xFtc555-1_CP003827.1_Chr8_QTL_test.vcf.gz
