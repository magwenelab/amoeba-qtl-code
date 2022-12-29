# Bring in needed mods
import numpy as np, pandas as pd
import scipy.stats as ss

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

def crypto_kruskal(site,pheno):
    refxpheno = np.array(pheno,dtype=float)[np.array(site)==0]
    altxpheno = np.array(pheno,dtype=float)[np.array(site)==1]
    return -np.log10(ss.kruskal(refxpheno,altxpheno)[1])

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
loci_path = "../../GENOTYPE/Bt22xFtc555-1_loci_cor.csv.gz"
chrommap_path = '../../GENOTYPE/H99_chrommap.csv'

## Bring in genotype data gather the segregants
## and show the frist five rows
geno = pd.read_csv(loci_path, index_col = 0)
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

## Set permutations, 
## Gather phenotype into an array
## make list of permutated data.
perms = 10000
pheno = Halo.loc[sorted(Map_pop),'Halo'].values
perms_pheno = [np.random.permutation(pheno) for i in range(perms)]

## Gather loci and intiate null distrbution
loci = geno[sorted(Map_pop)].drop_duplicates() 
null = []

## use a while loop to run permutaitons
i = 0
while len(null) < perms:
        
    ## Permute the phenotypic space
    null.append(loci.apply(allelic_manu,
                               args=[perms_pheno[i]],
                               axis=1).max())
    i += 1 ## Add one to i

## Save out dataframe
pd.DataFrame(null).to_csv('../../NOTES/Amoeba_halo_permutations.csv.gz',
                          index=False)
