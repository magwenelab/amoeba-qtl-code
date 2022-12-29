## Bring in needed mods
import pandas as pd, numpy as np, glob
import scipy.stats as ss, seaborn as sns
from matplotlib import pyplot as plt
#import matplotlib.gridspec as gridspec

## Set the color map to be used
qtlcolormap = ['#006BA4', '#FF800E', 
               '#ABABAB', '#595959', 
               '#5F9ED1', '#C85200', 
               '#898989', '#A2C8EC', 
               '#FFBC79', '#CFCFCF']

## Set the QTL ylabel that we use often
pvalylable = '$-\log_{10}\,(p-value)$'


## Define ftns for use in QTL mapping
def minmax_nrom(x):
    """Min max normalize values in the array X"""
    return (np.array(x) - np.min(x))/(np.max(x)-np.min(x))

## Write ftn for loading in genoetype data
def loadvariants(path1,path2,chroms):
    """
    Loads in genotype variant paths given wild card paths 
    (PATH1 and PATH2) and a list of chromosomes (CHROMS).
    """
    ## Gather allelic depth
    dfs = []

    ## Iterate thru chorms
    for i,c in enumerate(chroms):
    
        ## Gather paths pased on wilde cars
        info_path = glob.glob(path1%c)[0]
        dept_path = glob.glob(path2%c)[0]
    
        ## Load info and depth
        info = pd.read_csv(info_path,index_col=0)
        dp = pd.read_csv(dept_path,index_col=0)
    
        ## Concat dataframes
        depth = pd.concat([info,dp],axis=1)
    
        ## Append
        dfs.append(depth)
    
    ## Concat rows
    return pd.concat(dfs,axis=0).reset_index()

## Load in gff file given a path
def loadgff(path,
            sep='\t',
            header=None,
            comment='#',
            genecol = 'Gene',
            strandcol = 'Strand',
            parentcol='Parent',
            genesplit = '-t26',
            parentsplit = 'Parent=',
            idsplit = 'ID=',
            names = ["Seqid", "Source", "Type", "Start",
                     "End", "Score","Strand", "Phase", "Attribute"],
            dtype = ["str","str","str","int","int",
                     "str","str","str","str"]):
    """
    Loads in and returns a gene features file (GFF) given a specified PATH.
    The GFF is parsed using the seperated specified in SEP.
    A default header of NONE is set and COMMENTs are assumed to start with 
    an # symbol. 
    
    New columns add to the GFF by this function are named in
    GENECOL, STRANDCOL, and PARENTCOL which name the 
    genes, strand, and parent transcripts per gene, respectively. 
    
    The strings set in GENESPLIT, PARENTSPLIT, and IDSPLIT are used
    for splitting the attributes field within the GFF file. 
    
    NAMES: the set of column names of the GFF file.
    DTYPE: the data types set per column in the GFF.
    """
    gff = pd.read_csv(path,comment=comment,sep=sep,header=header,
                  names=names,dtype=dict(zip(names,dtype)))

    gff[strandcol] = gff[strandcol].replace(dict(zip(['-','+'],[-1,1])))
    gff[parentcol] = [a.split(parentsplit)[-1].split(';')[0].split(idsplit)[-1] 
                      for a in gff.Attribute]
    gff[genecol] = [a.split(genesplit)[0] for a in gff.Parent]
    return gff


## Mean confidnece intreval for arrays
def mean_confidence_interval(data, confidence=0.95, axis=0):
    """
    Calcualtes the confidence interval (CI) give critical value in 
    CONFIDENCE around the mean of DATA along the axis in AXIS (defalut = 0). 

    Returns the mean, lower, and upper bounds of a 95 CI of DATA.
    """
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a,axis=axis), ss.sem(a,axis=axis)
    h = se * ss.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


## Make a chromosome map
def chrommap(df,
             chrom='Chrom',
             pos='Pos',
             names = ['Chrom','Length','Cumlen','Midpts']):
    """
    Generates a chromsome map for analysis and plotting.
    
    From the collection of variants in the dataframe DF, 
    this function generates a dataframe with the following columns (left to right) 
    the chromosomes number or name in CHROM, the length of chromsome, 
    the cumulative length, and the cumlative midpoints of the chromosome.
    
    Inputs:
    DF: a collection of genetic variants (rows) with assoicated columns CHROM and POS.
    CHROM: the name of a the chromosome column in the dataframe DF.
    POS: the name of the positonal column in the dataframe DF.
    NAMES: the list of columns for the chrommap (MAP) to be returned.
    
    Output:
    MAP: a dataframe containing chromsome lengths, cumulative lenghts,
         and cumulative midpts across called genetic variants.
    """
    if type(df)!=pd.core.frame.DataFrame:
        print('Input DF must be a dataframe.')
        return
    else:
        pass
    
    Map = df[[chrom,pos]].groupby([chrom]).agg('max').reset_index()
    
    if chrom!=names[0]:
        names[0] = chrom
    else:
        pass
    
    Map.columns = names[:2]
    Map[names[2]] = [0] + list(np.cumsum(Map.Length.values)[:-1])
    Map[names[3]] = (Map.Length.values/2) + Map.Cumlen.values
    return Map

## An allelic mann-u test
def allelic_manu(geno,pheno,x=0,y=1): ## geno and pheno need to be in the same position
    """
    Conducts a Mann-Whitney U test on the phenotype data in PHENO
    by the genotypic states listed in GENO.
    Accepted data types include pandas.Series or numpy.array.
    
    Assumes the order of data within GENO and PHENO are paired.
    Defaluts for the biallelic state in GENO are 0 and 1; set in X and Y.
    
    Returns the -log10 of the calculated p-value from the test.
    See ?scipy.stats.mannwhitneyu for further help.
    """
    pheno = np.array(pheno) ## sets the type for the data as an array
    geno = np.array(geno)
    ## Gather phenotypes by genotypes 
    ## Parse the genotype data as True for 0 and then 1 and 
    ## take the asscoiated index within the phenotypic data array
    ## Return p-value
    return -np.log10(ss.mannwhitneyu(pheno[(geno==x)],pheno[(geno==y)])[1])

## An allelic kruskal W test.
def allelic_kruskal(site,pheno,x=0,y=1):
    """
    Conducts a Kruskal-Wallis H-test on the phenotype data in PHENO
    grouped by the genotypic states listed in GENO.
    Accepted data types include pandas.Series or numpy.array.
    
    Assumes the order of data within GENO and PHENO are paired.
    Defaluts for the biallelic state in GENO are 0 and 1; set in X and Y.
    
    Returns the -log10 of the calculated p-value from the test.
    See ?scipy.stats.kruskal for additional help and information.
    """
    refxpheno = np.array(pheno,dtype=float)[np.array(site)==x]
    altxpheno = np.array(pheno,dtype=float)[np.array(site)==y]
    return -np.log10(ss.kruskal(refxpheno,altxpheno)[1])

## Linear regression of phenotype on genotype.
def allelic_ANOVA(site, pheno):
    """
    This regression is equivalent to one-way ANOVA with 2 groups.
    Regresses the phenotype in PHENO on the genetic variant coded in SITE.
    Assumes genotypes in SITE are biallelic and coded as 0 and 1.
    
    Returns a coefficient of determination (R) and F-statistic (F).
    """
    coding = np.array(site, np.float)
    pheno = np.array(pheno, np.float)
    
    meany = np.mean(pheno)
    meandummy = np.mean(coding)
    ctry = pheno - meany
    ctrdummy = coding - meandummy
    
    # regression coefficient and intercept
    b = np.dot(ctry, ctrdummy)/np.dot(ctrdummy, ctrdummy)
    intercept = meany - b * meandummy
    
    yhat = b * ctrdummy
    len_yhat = np.sqrt(np.dot(yhat,yhat))
    len_y = np.sqrt(np.dot(ctry,ctry))
    df_yhat = 1
    
    error = ctry  - yhat
    len_error = np.sqrt(np.dot(error,error))
    if abs(len_error**2) < 1e-5:
        raise Exception("Zero length error in ANOVA")
    df_error = len(pheno) - 2
    
    # coefficient of determination is R**2
    R = (len_yhat/len_y)**2
    # F-statistic
    F = (len_yhat**2/df_yhat) / (len_error**2/df_error)
    return F,R

## Calculate p-value
def association_logPval(site, pheno):
    """
    Given a bi-allelic genetic variant (SITE) and phenotype (PHENO)
    Calculates an F-statistic (see allelic_ANOVA) and returns 
    the associated -log10(p-value) from the linear regression of PHENO on SITE.
    """
    F,R = allelic_ANOVA(site, pheno)
    logP = np.log10(ss.f.sf(F, 1, len(pheno)-2))
    return -logP
 
## Conducts QTL mapping    
def QTLmap(df,ftn,mappop,pheno):
    """
    Given a dataframe (DF) - with rows representing variant sites
    and columns as sample genotypes (coded as 0 or 1) - cunducts 
    QTL mapping using the function specified in FTN to associated 
    bi-allelic genotypes of the mapping population (MAPPOP) with 
    their phenotypes in PHENO. This function assumes the array in 
    PHENO is ordered with samples in MAPPOP.
    """
    loci = df[mappop].drop_duplicates() ## Reduces the number of unique test sites
    loci['QTL'] = loci.apply(ftn,args=[pheno],axis=1) ## QTL map
    return df.merge(loci.T.drop_duplicates().T) ## return merged results

## subset a chromosome variants
def findchrom(df,c, chrom='Chrom', pos='Pos', bounds = None):
    """
    Splits a dataframe (DF) on the column CHROM where the 
    value of CHROM is C. DF must be a dataframe were rows 
    are genetic variant sites with columns named CHROM and POS 
    representing the chormosome and position of variants.
    
    BOUNDS (defalt = None) is a tuple of positional
    coordinates to restrict the returned variants of a chromsome. 
    By defalut the left and right bounds are set to zero
    and the maximum genetic position along a given chromosome in DF.
    """
    if bounds is None:
        bounds = (0,df[pos].max())
    
    return df[(df[chrom]==c) & 
              (df[pos]>=np.min(bounds)) & 
              (df[pos]<=np.max(bounds))].sort_values(pos)

## Generate a map; e.g. figure 1 from Roth et al. 2018
def haplotypemap(geno,mappop,
                 c=None,
                 plotseg=None,
                 ns=10,
                 ws=10000,
                 slide=2500,
                 savepath=None,
                 chrom='Chrom',
                 pos='Pos',
                 parent_labels = None,
                 bicolor=['tab:blue','tab:orange'],
                 marker ='.',
                 hspace = 0.4,
                 fs = 12,
                 figsize=(10,5),
                 cobuff = 0.05,
                 mydpi=200,
                 samplelabelsize=6,
                 y_jitter=0.05,
                 legendbox = (-0.07,0.5),
                 legendscale=6,
                 legendtitle='Genotype',
                 centromerepath = None):
    """
    Generates a haplotype map of variants in GENO displaying 
    a moving histogram of variants along a chromosome, 
    the haplotypes of segregants listed in PLOTSEG, 
    and a plot of recombination counts. 
    
    The moving histogram is calculated using a window size of WS, 
    shifting every SLIDE, left to right along a chromosome.
    
    C: the chromosome that will be displayed. 
       Default is None and a chromosome from GENO is randomly chosen
       
    PLOTSEG: A list of segregants whose haplotypes to plot.
             Default is NONE and NS are selegted randomly.
             
    SAVEPATH: Path to save generated figure
    CHROM: Name of chromosome column in GENO, default = 'Chrom'
    POS: Name of positonal column in GENO, default = 'Pos' 
    BICOLOR: Colors used to represent genotypes in GENO
    MARKER: Style of marker used in plot, default = '.'
    HSPACE: Space between figures, default = 0.2
    FS: Fontsize used in plotting
    FIGSIZE: Length and width of figure
    COBUFF: Y-axis limit buff for the crossover plot
    MYDPI: The dpi to save the figure
    SAMPLELABELSIZE: Fontsize to set sample labels
    
    Returns or saves a figure of the haplotype map.
    """
    ## Call set up figure instance, grid, and axes
    fig = plt.figure(figsize=figsize)
    fig.set_facecolor('w')
    gs = fig.add_gridspec(5,1)
    vax = fig.add_subplot(gs[0,:])
    aax = fig.add_subplot(gs[1,:])
    hax = fig.add_subplot(gs[2:-1,:])
    cax = fig.add_subplot(gs[-1,:])

    ## Gather chromosome genotypes
    if (c is None):
        c = np.min(np.random.choice(geno[chrom].unique(),1))
    else:
        pass
    ctemp = findchrom(geno,c,chrom=chrom,pos=pos)
    
    ## Plot histogram of snps
    plt.sca(vax)
    
    ## Generate windows for histogram
    wss = np.arange(ctemp[pos].min()-slide,ctemp[pos].max()+slide,slide)
    wx = []
    wy = []
    for w in wss:
        wc = ctemp[(ctemp[pos]>=w) & (ctemp[pos]<=w+ws)]
        wx.append(wc[pos].mean())
        wy.append(wc.shape[0])
    plt.fill_between(wx,wy,color='tab:green',alpha=0.2)
    
    xp,xl = plt.xticks()
    plt.xticks(xp,labels=[])
    plt.xlim(ctemp[pos].min(),ctemp[pos].max())
    plt.ylabel('# of variants\nper %s kb'%int(ws/1000),
               rotation=0,
               fontsize=fs,
               verticalalignment='center',
               horizontalalignment='right')
    
    ## plot chromosome allele frequencies
    plt.sca(aax)
    af = ctemp[mappop].mean(axis=1)
    plt.plot(ctemp[pos],af.values,marker,color='k',alpha=0.5,markersize=1)
    plt.ylim(0,1)

    ## Plot centromere if path to bed file is non None
    if centromerepath is not None:
        centromeres = pd.read_csv(centromerepath)
        x = np.arange(*findchrom(centromeres,c,pos='Start')[['Start','End']].values[0])
        y = np.ones(len(x))
        plt.fill_between(x,y,color='grey',alpha=0.3)


    plt.ylabel('Allele\nfrequency',rotation=0,
               fontsize=fs,
               verticalalignment='center',
               horizontalalignment='right')
    xp,xl = plt.xticks()
    plt.xticks(xp,labels=[])
    plt.xlim(ctemp[pos].min(),ctemp[pos].max())
    
    ## Gather segregants to plot
    if (plotseg is None) and (type(plotseg)!= list):
        plotseg = list(np.random.choice(mappop,ns,replace=False))
    
    ## Plot sample haplotypes
    plt.sca(hax)
    
    if parent_labels is None:
        parent_labels = sorted(np.unique(np.concatenate(ctemp[plotseg].values)))
    
    for i,s in enumerate(plotseg):
        stemp = ctemp[[pos,s]]
        for j,g in enumerate([0,1]):
            hap = stemp[(stemp[s]==g)]
            plt.plot(hap[pos],
                     np.ones(len(hap[s]))*i+np.random.normal(0,y_jitter,len(hap[s])),
                     marker,markersize=2,
                     color=bicolor[j],alpha=0.2,
                     label = parent_labels[j] if i == 0 else None)
    
    if (legendbox is not None):
        plt.legend(bbox_to_anchor=legendbox,
               ncol=1,markerscale=legendscale,
               title=legendtitle)
            
    xp,xl = plt.xticks()
    plt.xticks(xp,labels=[])
    plt.yticks(np.arange(0,len(plotseg),1),plotseg,fontsize=samplelabelsize)
    plt.ylabel('Example\nsegregant\nhaplotypes',rotation=0,
               fontsize=fs,
               verticalalignment='bottom',
               horizontalalignment='right')
    plt.xlim(ctemp[pos].min(),ctemp[pos].max())

    ## Plot crossover plot
    plt.sca(cax)
    cof = ctemp[mappop].diff(axis=0).abs().sum(axis=1)/len(mappop)
    plt.plot(ctemp[pos].values,cof,color='k',alpha=0.5)
    plt.ylim(-cobuff,np.max(cof)+cobuff)
    xp,xl = plt.xticks()
    plt.xticks(xp,labels=[int(x/1000) for x in xp],fontsize=fs)
    plt.xlabel('Chromosome %s coordinates (kb)'%c,fontsize=fs)
    plt.ylabel('Crossover\nfrequency',rotation=0,
               fontsize=fs,
               verticalalignment='center',
               horizontalalignment='right')
    plt.xlim(ctemp[pos].min(),ctemp[pos].max())
    plt.subplots_adjust(hspace=hspace)
    if savepath is not None:
        plt.savefig(savepath,dpi=mydpi,bbox_inches='tight')
    return fig

## find the peak of QTL
def findpeak(pvaldf, c=None, pval='QTL', pos='Pos', chrom='Chrom'):
    """
    Given a dataframe (PVALDF) which represents the results of 
    QTL mapping between a chromosomes genotypes and sample 
    phenotypes returns the peak positions of a QTL.
    
    PVALDF must have two columns whoes values represent the strength 
    in association (PVAL, default = QTL) and postion of a variant (POS).
    
    If the value of C is not NONE the dataframe, PVALDF will be split 
    by the values in the column CHROM.
    
    This function returns two coordiantes, representing the left and right
    most positons (in the column POS) of genetic variants in the peak 
    of the QTL -- defined by values that have the same max 
    strength in association --, the max strength in association, 
    and the dataframe index.
    """
    if (c is not None):
        df = findchrom(pvaldf,c,chrom=chrom)
    else:
        df = pvaldf
    dfmax = df[(df[pval] == df[pval].max())]
    return dfmax[pos].min(),dfmax[pos].max(),dfmax[pval].max(),dfmax[pval].idxmax()
 
## Emperically calculate a confidence interval (CI) for a QTL  
def bootci(df, ftn, mappop, pheno,
           nboots=500,
           alpha=5,
           Pos = 'Pos',
           bootnames = ['l','r','m','QTL']):
    """
    Given chromosome genetic varinats in DF - where rows are 
    variants and columns are sample genotypes (coded as 0 or 1) -
    calculate a CI (using the number of bootstraps set in NBOOTS) 
    around a QTL (using the mapping function, FTN) for the phenotype, 
    PHENO, given the level of significance in ALPHA.
    
    Assumes the sample names and phenotypes of MAPPOP and PHENO 
    are paired.
    
    Returns the left and right positions of the CI and the dataframe
    of the bootstraps with column names set in BOOTNAMES.
    """
    sampdf = pd.DataFrame(np.array(pheno,dtype=float),
                          columns=['Phenotype'],index=mappop)
    
    boots = []
    while len(boots)<nboots:
        randf = sampdf.sample(n=len(mappop),replace=True)
        ransegs = randf.index.tolist()
        ranpheno = randf.Phenotype.values
        
        qtl = QTLmap(df,ftn,ransegs,ranpheno)
        l,r,pval,ix = findpeak(qtl)
        
        m = df[(df[Pos]>=l) & (df[Pos]<=r)][Pos].median()
        
        boots.append((l,r,m,pval))
    
    boots = pd.DataFrame(boots,columns=bootnames)
    cil = np.percentile(boots.m.values,alpha/2)
    cir = np.percentile(boots.m.values,100-(alpha/2))
    return cil,cir,boots


## Set a ftn for making a permutation dataframe
def permute_progeny(progeny, n = 1000):
    """
    Permutes the segregants listed in PROGENY, N times.
    Returns a pandas dataframe of N rows with K columns, 
    where K = len(progeny), that contains the permuted list of progeny.
    """
    return pd.DataFrame([np.random.permutation(progeny) for i in range(n)])

    
## Establish a null distribution using permutation tests
def permute_QTL(df, ftn, permdf, pheno, alpha=5):
    """
    Given a dataframe (DF) - whoes rows are variants and columns are 
    sample genotypes (coded as 0 or 1) - cunduct permutation (with N = 1000) 
    tests to establish a null distribution of association for QTL 
    mapping experiments using the function in FTN to associated 
    bi-allelic genotypes of the mapping population (MAPPOP) 
    with their phenotypes in PHENO.
    
    DF: A pandas dataframe with rows as variant positions and columns as 
        mapping population genotypes. Coded as zero or one.
    
    FTN: Name of function for use in QTL mapping.
    
    PERMDF: A dataframe of permuted mapping population samples names.
            This datarame is N x K matrix with N is the number of permutations
            and K is the number of unique progeny names.
            Make this dataframe with the function, PERMUTE_PROGENY.
            
    PHENO: An array of phenotypes.
    
    ALPHA: The critical value for use in establishing threshold (Default = 5). 
    
    Returns the threshold value from permutaion tests given ALPHA and 
    the null distribution of assoication as a list stored in PERMS. 
    """
    pheno = np.array(pheno,dtype=float)
    perms = -np.ones(permdf.shape[0])
    
    assert permdf.shape[1] == len(pheno), "The number of phenotypes != number of samples"
    
    for i,j in permdf.iterrows():
        mappop = j.values
        perms[i] = QTLmap(df,ftn,mappop,pheno)['QTL'].max()
        
    return np.percentile(perms, 100 - alpha),perms
        
## Display a manhattan plot
def manhattan(df,
              pval='QTL',
              pos='Pos',
              chrom='Chrom',
              chrmap = None, 
              ax = None, 
              threshold = None,
              figsize=(12,3), 
              savepath = None,
              ylims=None,
              xlims=None,
              ms =1, 
              fs =12, 
              rast = True, 
              mydpi = 200, 
              ylabel = pvalylable,
              xlabel = 'Chromosome',
              cmap = None,
              chromlabels = None,
              linestyle=' ',
              marker='.'):
    """
    Generates a manhattan plot given a dataframe (DF) with a column of 
    -log10(p-values) (PVAL, default = QTL) per genetic site (rows) and 
    additional columns represeting the position (POS) 
    and chromosome (CHROM) of a variant. 
    
    CHROMMAP: If provided, is a dataframe of lengths, cumaltive lengths, and 
              cumulative midpts of per chromosome. If NONE, CHROMMAP will be 
              auto generated from DF, see ?QTLfun.chrommap
    AX: a figure axes for plotting the manhattan plot.
        If NONE, ax will be generated
    THRESHOLD: A value for the height of a horizontal line, representing the 
               signifance threshold of a QTL mapping experiment. 
               The line will be plotted if thrdhold is set to not FALES.
    FIGSIZE: pair of integers representing the length and width of the figure.
    SAVEPATH: path to save figure if not none.
    MS: markersize for plotting
    FS: fontsize
    RAST: Boolean variable setting rasterization of plot.
    MYDPI: dpi to save figure
    YLABEL and XLABEL: labels of the y and x axis of the manhattan plot
    
    Returns a plotting handle with a manhattan plot generated per CHROM 
    with the strength of association in PVAL.
    """
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=figsize)
        fig.set_facecolor('w')
        
    if chrmap is None:
        chrmap = chrommap(df,chrom=chrom)
    
    if cmap is None:
        cmap = qtlcolormap
        
    if chromlabels is None:
        chromlabels = chrmap[chrom].values
    
    plt.sca(ax) ## Set the plotting axis
    ## Make manhattan plot,itterating over chromosomes
    for i,(ci,c) in enumerate(chrmap.groupby(chrom)):
        temp = df[(df[chrom]==ci)]
        cumpos = c.Cumlen.min()
        plt.plot(temp[pos].values+cumpos,temp[pval].values,
                 marker=marker,linestyle=linestyle,
                 markersize=ms,rasterized=rast,
                 color=cmap[int(i%len(cmap))])
        
    if (threshold is not None):
        plt.hlines(threshold,0,cumpos+temp[pos].max(),color='k',
                   linestyle='--',alpha=0.25,linewidth=1)
        
    ## Adjust xticks and set lables
    plt.xticks(chrmap.Midpts.values,chromlabels,fontsize=fs)
    plt.xlabel(xlabel,fontsize=fs)
    
    plt.yticks(fontsize=fs)
    plt.ylabel(ylabel,fontsize=fs)
    
    if (ylims is not None):
        plt.ylim(ylims[0],ylims[1])
    
    if (xlims is not None):
        plt.xlim(xlims[0],xlims[1])
    
    if (savepath is not None):
        plt.savefig(savepath,dpi=mydpi,bbox_inches='tight')
    return ax

## Construct a genotype by phenotype dataframe
def qtldf(df,pheno,mappop,qtl='QTL',qtlix=None,
          names = ['Strain','Phenotype','Genotype']):
    """
    Create a dataframe of the mapping populations (MAPPOP) 
    phenotypes (PHENO) and genotypes at the peak of a QTL in DF.
    
    DF is a dataframe with rows representing genetic variants and columns 
    are sample (stored in MAPPOP) genotypes (coded as 0 or 1).
    
    The strength in association is stored in a column in DF named QTL.
    
    Returns an m x 2 dataframe, with an index of sample names
    where rows (m) represent samples values of the two columns which
    represent phenotypic and genotypic values. 
    """
    if qtlix is None:
        qtlix = df[qtl].idxmax()
    gxp = pd.DataFrame([mappop,
                        np.array(pheno,dtype=float),
                        np.array(df.loc[qtlix,mappop].values,dtype=float)],
                        index = names).T
    gxp.index = gxp[names[0]] ## Set index to sample names
    gxp.drop(names[0],axis=1,inplace=True)
    for c in gxp.columns:
        gxp[c] = gxp[c].apply(float)
    return gxp

## Plot genotypes by phenotypes at a QTL, as seen in Roth et al 2021.
def gpplot(gxp,Phenotype='Phenotype',Genotype = 'Genotype',Hue=None,
           parents = None,plabel=None,
           ylabel=None,xlabel="Parental Allele",annotation = True,
           figsize=(8,8),ax = None,fs=12,
           markersize=4,markersizemod=6,
           pcolors=['tab:blue','tab:orange']):
    """
    Generates a genotype by phenotype swarm plot from the genotype by phenotype
    dataframe, GXP. This dataframe has two columns, PHENOTYPE and GENOTYPE,
    which represent the phenotypes and genotypes of samples at a QTL. 
    The index of this dataframe are the sample (or progeny) names.
    
    PARENTS: a sorted list of the two parent strain names used in the cross.
             This function assumes this list matches allelic values of 0 and 1.
    PLABEL: a sorted list of the two parent strain names used in legned
    
    YLABEL and XLABEL: The labels of the y and x axis (respectively), 
                       defaluts are None and "Parental Allele".
    ANNOTATION: A boolean variable to annotate plot with the 
                coefficnet of determiantion.
    FIGSIZE: the size (length and width) of a figure.
    AX: figure axis for plotting.
    FS: The fontsize of lables and annotation.
    MARKERSIZE: The markersize of used in plotting.
    MARKERSIZEMOD: The additional size to add to markersize.
    PCOLORS: Colors of parental alleles.
    """
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=figsize)
        fig.set_facecolor('w')
        
    plt.sca(ax)
        
    sns.swarmplot(x = Genotype, y = Phenotype,
                  data = gxp,
                  color = "grey",
                  size=markersize,
                  alpha=0.5,hue=Hue);
    
    sns.regplot(x = Genotype, y = Phenotype,data = gxp,
                ci=False,
                line_kws = {"alpha":0.8, "color":"red"}, 
                scatter_kws = {"alpha":0});

    plt.xticks([0,1], plabel, fontsize = fs);
    plt.xlabel(xlabel, fontsize = fs);
    
    plt.yticks(fontsize = fs);
    plt.ylabel(ylabel,fontsize=fs)

    plt.plot(gxp.groupby(Genotype).mean().index,
             gxp.groupby(Genotype).mean()[Phenotype],
             'k_',markersize=markersize+markersizemod)
    
    if (parents is not None) and (plabel is not None):
        [plt.plot(gxp.loc[p,Genotype],gxp.loc[p,Phenotype],'^',
              color = pcolors[i],alpha=0.8,label = plabel[i],
              markersize=markersize+markersizemod) 
              for i,p in enumerate(parents)]
        
    if annotation:
        F,R = allelic_ANOVA(gxp[Genotype].values,gxp[Phenotype].values)
        
        max90 = gxp[Phenotype].max()*(0.9)
        plt.annotate('$R^2$ = %s'%np.round(R,3),(0.5,max90),
                    va='center',ha='center',fontsize=fs)
        #plt.title(x=0.5,y=0.85,
        #      label='$R^2$ = %s'%np.round(R,3),
        #      va='center',ha='center',fontsize=fs)
    
    return ax

## Temporal QTL mapping
def longitudinal_regression(geno,pheno,dr=1,de=2):
    """
    Regress a temporal (or longitudinal) phenoytpe onto a genotype.
    
    GENO: An 1 x m pandas series representing a singluar, bi-allelic 
          genetic site coded with 0 and 1 (or -1 and 1) 
          and m are sample names.
          
    PHENO: An m x p dataframe where m are samples -- ordered
           top to bottom to match their apperance in GENO, 
           left to right -- and p are ordered, continues 
           logitudinal index of the phenotypes from left to right.
    
    Returns the coefficent of determination (RSQ) and the 
    -log10(p-value) (LOGP) across the longitudinal axis.
    """
    betadf = pd.DataFrame(
                [np.ones(len(geno)),geno.values],
                   columns=geno.index,
                index=['b0','b1']).T
        
    res = np.linalg.lstsq(betadf.values, pheno,rcond=-1)
    coeff = res[0]
    yhat = coeff[0] + np.outer(betadf['b1'].values,coeff[1])
    
    sse = np.sum((pheno.values-yhat)**2,axis=0)
    ssy = np.sum((pheno - pheno.mean())**2,axis=0).values
    rsq = (ssy - sse)/ssy
    
    mse = sse/(len(geno)-de)
    msr = (ssy-sse)/dr
    fratio = msr/mse
    logP = -np.log10(ss.f.sf(fratio, dr, len(geno)-de))
    
    return rsq,logP

## Wrapper of temporal QTL mapping
def longitudinal_qtl(geno,pheno,mappop,ftnix=None,dr=1,de=2):
    """
    Conduct a logitudinal regression to associate a functional
    phenotype with a genotype.
    
    GENO: An n x m dataframe where n are the representation of 
          genetics variants and m are sample names. 
          The values of GENO are assumed to be bi-allelic as
          either 0 and 1 or -1 and 1. 
    
    PHENO: An m x p dataframe were m are samples and p 
           are ordered phenotypes. The row index of PHENO must 
           be the sample names in MAPPOP.
    
    MAPPOP: A list of segregants/samples to use in mapping.
    
    FTNIX: An ordered list of index of the function valued phenotype.
           Default is NONE and taken from column names in PHENO.
    
    DR and DE: Degrees of freedom for the regressor and error vectors. 
    
    Returns two dataframes created from merger with the GENO dataframe
    that contain the longitudinal coefficient of determination (GRSQ)
    and the longitudinal strength of association (GLOG), respectively.
    """
    if ftnix is None:
        ftnix = pheno.columns
        
    temp = pheno.loc[mappop,ftnix]
    loci = geno[mappop].drop_duplicates()
    
    res = loci.apply(longitudinal_regression,args=[temp,dr,de],axis=1)
    
    rsq = pd.concat([loci,pd.DataFrame([a[0] for a in res.values],
                     index=loci.index,columns=ftnix)],axis=1)
    
    logp = pd.concat([loci,pd.DataFrame([a[1] for a in res.values],
                      index=loci.index,columns=ftnix)],axis=1)
    
    grsq = geno.merge(rsq)
    glog = geno.merge(logp)
    
    return grsq,glog

## Return the max and min acorss longitudinal dataframe
def longitudinal_stats(glog,ftnix,chrom='Chrom',pos='Pos'):
    """
    Given an (n x m) genotype dataframe GLOG -- with rows (n) 
    representing bi-allelic variant sties and columns (m) 
    rerpeseting samples, chromosomes (CHROM), and positon 
    of variants (POS) and the strength of association 
    (-log10(p-value)) across a functional axis (FTNIX) 
    returns the average, maximum, and final strength in 
    association across FTNIX per variant (n) appended to 
    GLOG as columns along with the chromosome (CHROM) and 
    positional (POS) coordinates columns per variant. 
    """
    glog['Average'] = glog[ftnix].mean(axis=1)
    glog['Maximum'] = glog[ftnix].max(axis=1)
    glog['Final'] = glog[ftnix[-1]]
    return glog[[chrom,pos]+['Average','Maximum','Final']].copy()

## Returns the max indexes along the chromosome positon and 
## Funcitonal index for QTL 
def findmaxix(logp,ftnix,
              chrom='Chrom',pos='Pos',
              names=['Chrom','Pos','Ftnix','QTL']):
    """
    Given a longitudinal dataframe of QTL mapping results (LOGP) with a 
    functional index of FTNIX (such as time points) and columns representing
    the chromosomes (CHROM) and variant position (POS) this funciton returns 
    a dataframe with both the index in FTNIX and POS that maximize the 
    strength in associaiton per chromosome in CHROM. The NAMES argument is a 
    list of strings for the retuned summary dataframe. 
    """
    res = []
    for c,ctemp in logp.groupby(chrom):
        posmix = ctemp[ftnix].max(axis=1).idxmax()
        ftnmix = ctemp[ftnix].max(axis=0).idxmax()
        maxp = ctemp[ftnix].max().max()
        res.append((c,ctemp.loc[posmix,pos],ftnmix,maxp))
    return pd.DataFrame(res,columns=names)

## manhattan heat map
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
def manhattan_heatmap(logp,ftnix,ax=None,chrmap=None,
                      chrom='Chrom',pos='Pos',
                      colormap= "YlGnBu",ftnsplit = '_',
                      levels = None,spines=True,
                      alpha=0.6,linemin=0,linealpha=0.5,
                      figsize=(8,3),fs=12,
                      ylabel=None,xlabel=None,
                      cbar=True,pvalylable=pvalylable,
                      cbarw = 5,cbarh = 100,fc ='w'):
    """
    Creates a Manhattan heatmap for longituidnal QTL analysis.
    """
    ## Make a fig axis if none is fed into the ftn
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=figsize)
        fig.patch.set_facecolor(fc)
        
    ## Convert the functional index into an array
    if type(ftnix[0]) == str:
        ftnx = np.array([float(x.split(ftnsplit)[-1]) for x in ftnix])
    else:
        ftnx = np.array(ftnix)
        
    ## Make sure ftnix is a list
    if type(ftnix) != list:
        ftnix = list(ftnix)
    
    ## Set pvalue levels
    if levels is None:
        levels = np.arange(np.round(logp[ftnix].max().max())+1)

    ## Make a 
    if chrmap is None:
        chrmap = chrommap(logp)
        
    chrlist = chrmap.Chrom
    cumlen = chrmap.Cumlen
    clens =  chrmap.Length
    midpts = chrmap.Midpts
    
    ## Iterate over chromoosmes
    for c,chrom in enumerate(chrlist):
        cumpos = cumlen[c]
        temp = logp[(logp.Chrom==chrom)][ftnix+[pos]]
        chpos = temp.Pos.values + cumpos
        
        ## Plot heat map
        cs = ax.contourf(chpos,ftnx,
                         temp[ftnix].T.values,
                         cmap=colormap,
                         levels=levels,
                         extend='max',
                         alpha=alpha);
        
    ## Set xlim    
    plt.xlim(0,temp.Pos.max() + cumpos)
    
    plt.sca(ax)
    ## plot vertical lines seperating chromosomes
    plt.vlines(np.array(clens) + np.array(cumlen),
               linemin,np.max(ftnx),
               color='k',linewidth=1,
               linestyle='--',alpha=linealpha);
    
    ## Set xticks
    plt.xlabel(xlabel,fontsize=fs)
    plt.xticks(midpts,labels = chrlist, fontsize=fs)
    plt.ylabel(ylabel,fontsize=fs)
    plt.yticks(fontsize=fs)
    
    ## Turn off spines if true
    if spines:
        [ax.spines[k].set_visible(False) 
         for k in ['top','right']];
    
    ## add the color bar    
    if cbar:
        axins = inset_axes(ax,
           width="%s"%(cbarw)+'%', # width = % of parent_bbox width
           height="%s"%(cbarh)+'%', # height = % of parent_bbox height
           loc=6,
           bbox_to_anchor=(1.05, 0., 1, 1),
           bbox_transform=ax.transAxes,
           borderpad=0,)

        cbarp = plt.colorbar(cs, cax=axins)
        plt.ylabel(pvalylable,fontsize=fs)
        cbarp.set_ticks(levels[::2])

    return ax


def modelprint(res,png,
               x = 0.01,y=0.05,mydpi=80,
               figsize=(8,5),fs=10,fontprop='monospace'):
    """
    Prints an OLS model result summary page (RES) to figure for saving as PNG.
    """
    plt.rc('figure', figsize=figsize, facecolor='w');
    plt.text(x, y, str(res.summary()), 
         {'fontsize': fs}, fontproperties = fontprop);
    plt.axis('off');
    plt.savefig(png,bbox_inches='tight',dpi=mydpi);
    plt.close()
    
    
## Ftns for LDA analysis    
def return_mean_vector(X, y):
    """
    Returns the mean vector of phenotypes in X based on unique, sorted classes in Y.
    """
    return np.array([np.mean(X[y==c],axis=0) for c in np.unique(y)])


def within_scatter(X, y):
    """
    For LDA on X using unique classes in Y, calculates the within scatter matrix S_W.
    """
    class_labels = np.unique(y)
    n_classes = class_labels.shape[0]
    n_features = X.shape[1]
    mean_vectors = return_mean_vector(X, y)
    S_W = np.zeros((n_features, n_features))
    for cl, mv in zip(class_labels, mean_vectors):
        class_sc_mat = np.zeros((n_features, n_features))                 
        for row in X[y == cl]:
            row, mv = row.reshape(n_features, 1), mv.reshape(n_features, 1)
            class_sc_mat += (row-mv).dot((row-mv).T)
        S_W += class_sc_mat                           
    return S_W


def between_scatter(X, y):
    """
    For LDA on X using unique classes in Y, calculates the between scatter matrix S_B.
    """
    overall_mean = np.mean(X, axis=0)
    n_features = X.shape[1]
    mean_vectors = return_mean_vector(X, y)    
    S_B = np.zeros((n_features, n_features))
    for i, mean_vec in enumerate(mean_vectors):  
        n = X[y==i+1,:].shape[0]
        mean_vec = mean_vec.reshape(n_features, 1)
        overall_mean = overall_mean.reshape(n_features, 1)
        S_B += n * (mean_vec - overall_mean).dot((mean_vec - overall_mean).T)
    return S_B


def LDAeigenvalues(X,y):
    """
    Performs LDA on X given descret classes in Y.
    Returns the calulcated eigenvalues (and vectors) calculated from within and between scatter matricies. 
    """
    return np.linalg.eig(np.linalg.inv(within_scatter(X, y)).dot(between_scatter(X, y)))


def LDAeigenpairs(X,y):
    """
    Returns the sorted eigenvalues and assocated eigenvectors from LDA on X given descret classes in Y.
    """
    eigenvals, eigenvecs = LDAeigenvalues(X,y)
    eigenpair = [(np.abs(eigenvals[i]), eigenvecs[:,i]) for i in range(len(eigenvals))]
    eigenpair = sorted(eigenpair, key=lambda k: k[0], reverse=True)
    return np.array([k[0] for k in eigenpair]),np.array([k[1] for k in eigenpair])
    

def LDAcomponents(X, y, n=2):
    """
    Returns the N stacked eigenvectors from LDA on X given the descret classes in Y.
    Eigenvectors are sorted and returned in order of their maximum real eigenvalue.
    """
    eigenvals, eigenvecs = LDAeigenpairs(X,y)
    return np.hstack([eigenvecs[i].reshape(len(eigenvals), 1) for i in range(n)])


def variance_explained(eigenvals):
    """
    Return a np.array of summed ratios of eigan values in EIGVALS.
    Only returns sum of the real parts of each eigan value.
    """
    return np.real(np.array(eigenvals))/np.sum(np.real(np.array(eigenvals)))


def LDAmax(y,X):
    """
    Returns the maximum eigenvalue from LDA on data in X given descret class in Y.
    """
    eigenvals, eigenvecs = LDAeigenpairs(X,np.array(y).astype(int))
    return np.max(eigenvals)
    
    