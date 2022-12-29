#!/usr/bin/env python
# coding: utf-8

# In[16]:

## Set strain
si =1 
parental_pmy = ['PMY2649', 'PMY2650']
samplename = parental_pmy[si]


# In[1]:

## Bring in mods and write ftns for analysis
import pandas as pd, numpy as np, os, glob
from Bio import SeqIO

def dotfill(r,a,fill='.'):
    return ''.join([fill for i in range(abs(len(a)-len(r)))])

def scorealign(r,a,fill='.',verbose=False):
    score = 0
    ns = []
    for i in range(len(r))[1:-1]:
        temp = r[:i] + dotfill(r,a) + r[i:]
        c = 0
        for j,n in enumerate(a):
            if n == temp[j]:
                c = c + 1
        if c > score:
            ns = temp
            score = c
            if verbose:
                print(ns,c)
    if len(ns) == 0:
        #print('Error empty sequence\n%s\n%s'%(r,a))
        ns = r[0]+dotfill(r,a)+r[-1]
    return ns

def calcpos(seq,start,end,pad=0,fill='.'):
    c = 0;
    pos = np.arange(start,end+pad)
    poss= []
    if fill in seq:
        for n in seq:
            if n != fill:
                poss.append(pos[c])
                c = c + 1
            else:
                poss.append(pos[c])
    else:
        poss = pos
    return poss

def realign(ref,alt,verbose=False):
    if len(alt) > len(ref):
        newref = scorealign(ref,alt,verbose=verbose)
        newalt = alt
    elif len(ref) > len(alt):
        newalt = scorealign(alt,ref,verbose=verbose)
        newref = ref
    else:
        assert len(alt) == len(ref)
        newref = ref
        newalt = alt
    return newref,newalt

rdict = dict(zip(['three_prime_UTR',
                  'five_prime_UTR',
                  'CDS'],
                 [3,5,0]))

def gvtodf(gv,refseqdf,gffdf,rs='REF',fill='.',rdict=rdict,verbose=False):

    refs = []
    alts = []
    poss = []
    for i,j in gv.iterrows():
    
        ref = j.Ref
        alt = j.Alt.split('.')[int(j.Sample)]
        start = j.Pos
        end = start + len(ref)
    
        newref,newalt = realign(ref,alt)
        newpos = calcpos(newref,start,end)
    
        refs.append(newref)
        alts.append(newalt)
        poss.append(newpos)
    
    pos = np.concatenate(poss)
    ref = [n for n in ''.join(refs)]
    alt = [n for n in ''.join(alts)]

    df = pd.DataFrame([pos,ref,alt],
                      index=['Pos','Ref','Alt']).T
    df['Isvar'] = 1

    check = refseqdf[(refseqdf.Pos.isin(df.Pos))]
    failed = []
    for i,j in check.iterrows():
        
        tocheck = [s for s in 
                   df[(df.Pos==j.Pos)].Ref.tolist() 
                   if (s != fill)]
        if len(tocheck) > 1: ## Should be one base
            failed.append(j.Pos)
        elif not tocheck[0] == j[rs]:
            failed.append(j.Pos)
            
    #if (len(failed) > 0) & verbose:
        #print('Reference alignment failed')
    
    dfadd = refseqdf[~(refseqdf.Pos.isin(df.Pos))].copy()
    dfadd['Ref'] = dfadd[rs]
    dfadd['Alt'] = dfadd[rs]
    dfadd['Isvar'] = 0
    dfadd.drop(rs,axis=1,inplace=True)

    res = pd.concat([dfadd,df]).sort_values('Pos')
    res['Strand'] = gffdf.Strand.min()
    res['Type'] = -1
    res['Phase'] = -1

    for i,j in gffdf[(gffdf.Type!='gene')].iterrows():
    
        nb = np.arange(j.Start,j.End+1)
        res.loc[(res.Pos.isin(nb)),'Type'] = j.Type
    
        if j.Type == 'CDS':
            res.loc[(res.Pos.isin(nb)),'Phase'] = int(j.Phase)
    
    assert -1 not in res[(res.Type=='CDS')].Phase.tolist()
    
    res.Type.replace(rdict,inplace=True)

    return res,failed


# In[2]:


chrommap_path = '../../GENOTYPE/H99_chrommap.csv'
chrommap = pd.read_csv(chrommap_path)


# In[3]:


geno_path = '../../GENOTYPE/Bt22xFtc555-1_genotypes.csv.gz'
genos = pd.read_csv(geno_path,index_col=0)
genos = genos.merge(chrommap[['Chrom','Seqid']])


# In[5]:


if not 'Ref' in genos.columns:
    genos['Ref'] = [a.split('.')[0] for a in genos.Alleles]
    
if not 'Alt' in genos.columns:
    genos['Alt'] = genos['Alleles']


# In[6]:


## Bring in GFF file
gffpath = '../../REF/FungiDB-46_CneoformansH99.gff.gz'
names = ["Seqid", "Source", "Type", "Start", "End", "Score", 
         "Strand", "Phase", "Attribute"]
dtype = ["str","str","str","int","int","str","str","str","str"]

gff = pd.read_csv(gffpath,comment='#',
                  sep='\t',header=None,
                  names=names,dtype=dict(zip(names,dtype)))

gff['Strand'] = gff['Strand'].replace(dict(zip(['-','+'],[-1,1])))
gff['Parent'] = [a.split('Parent=')[-1].split(';')[0].split('ID=')[-1] 
                 for a in gff.Attribute ]

gff['Gene'] = [a.split('-t26')[0] for a in gff.Parent]
gff = gff.merge(chrommap)


# In[7]:


## Gather genes
genes = gff[(gff.Type=='gene')].sort_values(['Seqid','Start','Strand']).copy()

genes['Description'] = [a.split('description=')[-1].split('%2C')[0] 
                        for a in genes.Attribute]

# In[8]:


## Make list of descriptions and features
descriptions = ['hypothetical protein',
                'unspecified product',
                'conserved hypothetical protein']

foi = ['gene','three_prime_UTR','five_prime_UTR','CDS']


# In[9]:


## Bring in reference file
refpath = '../../REF/FungiDB-46_CneoformansH99_Genome.fasta'
REF = [s for s in SeqIO.parse(refpath,format='fasta')]


# In[17]:


savepath = '../../GENOTYPE/GENES/'
savepath_sample = savepath+'%s/'%samplename

if not os.path.exists(savepath):
    os.mkdir(savepath)

if not os.path.exists(savepath_sample):
    os.mkdir(savepath_sample)
    


# In[18]:


info_cols = [c for c in genos.columns if c[:3]!='PMY']


# In[20]:


## Is the GFF zero based?
zb = False

## Do you want to print everything 
## to the screne
## like a madman?!?!
verbos = True


# In[22]:


## RUN ALL THE GENES WITH CDS
#chroms = [1, 2, 3, 5, 6, 10, 11]
#run_parents = gff[(gff.Type=='CDS') & (gff.Chrom.isin(chroms))].Parent.unique()
run_parents = gff[(gff.Type=='CDS')].Parent.unique()

# In[23]:


novars = []
notype = []
sofail = []
for ssk1name in run_parents:
    genesave = savepath_sample+ssk1name+'.csv.gz'
    if os.path.exists(genesave):
       continue 
    ssk1df = gff[(gff.Type.isin(foi)) & 
                 (gff.Parent==ssk1name)
                ].sort_values('Start').copy()
    
    ssk1df['Start'] = ssk1df['Start'] + (1 if zb else 0)
    
    assert len(ssk1df.Strand.unique())==1
    if len(ssk1df.Type.unique())!=3:
        notype.append(ssk1name+'\n')
        continue

    sample = genos.loc[genos[samplename].dropna().index,
                       info_cols+[samplename]]
    
    sample['Sample'] = sample[samplename]
    ## Locate genetic variants within a gene
    gv = sample[(sample.Seqid==ssk1df.Seqid.min()) & 
            (sample.Pos>=ssk1df.Start.min()) & 
            (sample.Pos<=ssk1df.End.max())]
    if gv.shape[0] == 0:
        novars.append(ssk1name+'\n')
        continue
    ## From the referecne take the sequence of interest
    gene_seq = [s.seq[ssk1df.Start.min()-(0 if zb else 1):ssk1df.End.max()] 
            for s in REF if s.id == ssk1df.Seqid.min()][0]
    ## Make a dataframe from these sequences
    refseqdf = pd.DataFrame([np.arange(ssk1df.Start.min()+(1 if zb else 0),
                                       ssk1df.End.max()+1),
                     [n for n in str(gene_seq)]],
                            index=['Pos','REF']).T;
    ## Make gene variant dataframe
    res,failed  = gvtodf(gv,refseqdf,ssk1df,verbose=verbos)
    res['Gene'] = ssk1name
    res['Sample'] = samplename
    ## Save and print results
    res.drop_duplicates().to_csv(genesave,index=False)
    if verbos and (len(failed)>0):
        sofail.append(genesave)

# In[ ]:


failed_path = savepath_sample+'FAILED.csv'
if (len(sofail) > 0):
   open(failed_path,'w').writelines(sofail)
    
novars_path = savepath_sample+'NOVARS.csv'
if (len(novars) > 0):
   open(novars_path,'w').writelines(novars)

notype_path = savepath_sample+'NOTYPE.csv'
if (len(notype) > 0):
   open(notype_path,'w').writelines(notype)
