# Bring in needed mods
import numpy as np, pandas as pd
##from matplotlib import pyplot as plt
import scipy.stats as ss, seaborn as sns

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

## Define parental labels
parents = ['Bt22','Ftc555-1']

## Set paths to phenotype data
## 1) the parental data (parents)
## 2) the initial set of segregants (old)
## 3) the latest set of segregants (new)
## Note, these data where previously processed
## And the halo per replicates / experimetns where calcualted
parent_data_path = "../../PHENOTYPE/AMOEBA/HALOS/bt22xftc555_halo_new_parents.csv"
Old_data_path = "../../PHENOTYPE/AMOEBA/HALOS/bt22xftc555_halo_old.csv"
New_data_path = "../../PHENOTYPE/AMOEBA/HALOS/bt22xftc555_halo_new.csv"

## set path to genotype data and chromosome map
loci_path = "../../GENOTYPE/Bt22xFtc555-1_loci.csv.gz"
chrommap_path = '../../GENOTYPE/H99_chrommap.csv'

## Bring in genotype data gather the segregants
## and show the frist five rows
geno = pd.read_csv(loci_path, index_col = 0)
to_drop = geno[(geno.Chrom==8) & 
               (geno.Pos.isin([726391, 731769]))].index
geno.drop(to_drop,axis=0,inplace=True)
geno_segs = [s for s in geno.columns if s[:3] == "PMY"]

## Gather parental data
## and average by replicate per experiment
PNd = pd.read_csv(parent_data_path)
PNda = PNd.groupby(["PMY","Z"]).mean().reset_index().groupby("PMY").mean().reset_index()
PNda.drop("Z", axis = 1, inplace = True)

## Gather PMY numbers of paretns
Parent_PMY = PNda.PMY.tolist()

## Gather the latest data and average across repliacates, view head
Nd = pd.read_csv(New_data_path)
Nda = Nd[~(Nd.isin(Parent_PMY))].groupby("PMY").mean().reset_index()

## For the inital (old) data load in and average replicates
Od = pd.read_csv(Old_data_path)
Oda = Od[~(Od.isin(Parent_PMY))].groupby("PMY").mean().reset_index()

## Concatonate resluts
Halo = pd.concat([PNda, Oda, Nda]).reset_index(drop = True)

## Convert inches to cm
Halo["Halo_in"] = Halo.Halo
cm_con = 2.54
Halo["Halo"] = (np.sqrt(Halo.Halo_in.values)*cm_con)**2

## Set index
Halo.index = Halo.PMY

## Check the number of times each sample appears in dataframe
pmy, pn = np.unique(Halo.PMY.values, return_counts = True)

if np.max(pn)>1:
    print(pmy[(pn>1)])
    
## Gather segregants with both genotype and phenotype data
Map_pop = [s for s in Halo.PMY.unique() if s in geno_segs]
len(Map_pop) ## print # of segregants

## Set number of bootstraps, 
perms = 300

## Gather loci and intiate distrbution
chrom = geno[(geno.Chrom==8)][sorted(Map_pop)+['Pos']].copy().reset_index(drop=True)
bsp = []

## use a while loop to run bootstraps
i = 0
while len(bsp) < perms:
    
    bsegs = np.random.choice(Map_pop,len(Map_pop))
    bpheno = Halo.loc[bsegs,'Halo'].values
    bloci = chrom[bsegs].drop_duplicates()
    
    bloci['Pval'] = bloci.apply(allelic_manu,args=[bpheno],axis=1)
        
    bres = chrom.merge(bloci.T.drop_duplicates().T)
    bres.sort_values('Pval',inplace=True)
    
    bpos = bres[(bres.Pval==bres.Pval.max())].Pos
    bsp.append((bpos.median(),bpos.mean(),bpos.min(),bpos.max(),bres.Pval.max()))
    i += 1
    
bdf = pd.DataFrame(bsp,columns=['Median','Mean','Min','Max','Pval'])

## Save out dataframe
bdf.to_csv('../../NOTES/Amoeba_halo_bootstraps.csv.gz',
                          index=False)