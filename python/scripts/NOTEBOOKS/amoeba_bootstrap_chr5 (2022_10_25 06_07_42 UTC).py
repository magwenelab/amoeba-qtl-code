# Bring in needed mods
import numpy as np, pandas as pd, scipy.stats as ss

## Define ftns for use in QTL mapping
def allelic_manu(geno,pheno,x=0,y=1): ## geno and pheno need to be in the same position
    """
    Conducts a Mann-Whitney U test on the phenotype data in PHENO
    by the genotypic states listed in GENO.
    
    Assumes the order of data within GENO and PHENO are paired.
    Defaluts for the biallelic stat in GENO are 0 and 1; set in X and Y.
    
    Returns the -log10 of the calculated p-value
    """
    pheno = np.array(pheno) ## sets the type for the data as an array
    geno = np.array(geno)
    ## Gather phenotypes by genotypes 
    ## Parse the genotype data as True for 0 and then 1 and 
    ## take the asscoiated index within the phenotypic data array
    ## Return p-value
    return -np.log10(ss.mannwhitneyu(pheno[(geno==x)],pheno[(geno==y)])[1])

## Set the random seed
np.random.seed(71191)

## Set paths to phenotype data by genotype data and bring into environment
data_path = '../../NOTES/Bt22xFtc555-1_QTL8_genotype_phenotype.csv'
qtl = pd.read_csv(data_path,index_col=0)

## Bring in genotype data gather the segregants
## set path to genotype data and chromosome map
loci_path = "../../GENOTYPE/Bt22xFtc555-1_loci.csv.gz"
geno = pd.read_csv(loci_path, index_col = 0)

## Gather segregants and check length
geno_segs = [s for s in geno.columns if s[:3] == "PMY"]
assert len(geno_segs)>0

## Gather segregants with Ftc555-1 genotype
## at chromoosme 8 and check lenth of segregnats
Map_pop = [s for s in 
           qtl[(qtl.GT==1)].index.values if s in geno_segs]
assert len(Map_pop) > 0 
assert len(Map_pop) < qtl.dropna().shape[0]

## Set number of bootstraps
perms = 300

## Gather loci and intiate distrbution
chrom = geno[(geno.Chrom==5)][sorted(Map_pop)+['Pos']
                             ].copy().reset_index(drop=True)
bsp = []

## use a while loop to run bootstraps
i = 0
while len(bsp) < perms:
    
    bsegs = np.random.choice(Map_pop,len(Map_pop))
    bpheno = qtl.loc[bsegs,'Halo'].values
    bloci = chrom[bsegs].drop_duplicates()
    
    bloci['Pval'] = bloci.apply(allelic_manu,
                                args=[bpheno],
                                axis=1).replace(np.inf,np.nan)
        
    bres = chrom.merge(bloci.T.drop_duplicates().T)
    bres.dropna(axis=0,how='any',inplace=True)
    #bres.sort_values('Pval',inplace=True)
    
    bpos = bres[(bres.Pval==bres.Pval.max())].Pos
    bsp.append((bpos.median(),bpos.mean(),bpos.min(),bpos.max(),bres.Pval.max()))
    i += 1
    
bdf = pd.DataFrame(bsp,columns=['Median','Mean','Min','Max','Pval'])

## Save out dataframe
bdf.to_csv('../../NOTES/Amoeba_halo_bootstraps_chr5.csv.gz',
                          index=False)