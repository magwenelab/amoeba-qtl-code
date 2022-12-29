#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd, numpy as np, glob
from Bio.Seq import Seq

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def makeorf(df,fill='.'):
    ref = Seq(''.join(df[(df.Ref!='.')]['Ref'].tolist()))
    alt = Seq(''.join(df[(df.Alt!='.')]['Alt'].tolist()))
    
    if df.Strand.min() < 0:
        ref = ref.reverse_complement()
        alt = alt.reverse_complement()
        
    return ref,alt


# In[3]:


chrommap_path = '../../GENOTYPE/H99_chrommap.csv'
chrommap = pd.read_csv(chrommap_path)


# In[4]:


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

genes = gff[(gff.Type=='gene')].sort_values('Start').copy()

genes['Description'] = [a.split('description=')[-1].split('%2C')[0] 
                        for a in genes.Attribute]

#genes['Chromosome'] = [int(a[-2:]) for a in genes.Chrom]

genes.sort_values(['Chrom','Start'],inplace=True)


# In[5]:


novars_paths = glob.glob('../../GENOTYPE/GENES/PMY*/NOVARS.csv')

novars = np.unique(np.concatenate([pd.read_csv(p,header=None)[0].tolist() for p in novars_paths]))

no_vars = []
for nv in novars:
    temp = gff[(gff.Parent==nv) & (gff.Type=='CDS')]
    
    el = np.sum(temp['End'] - temp['Start']+1)/3
    refs = 1
    alts = 1
    nonsyn = 0
    nvars = 0
    utr5 = 0
    utr3 = 0
    exon = 0
    intron = 0
    chrom = temp.Chrom.min()
    descr = genes[(genes.Gene == temp.Gene.min())].Description.min()
    gstart = genes[(genes.Gene==temp.Gene.min())].Start.min()
    gend = genes[(genes.Gene==temp.Gene.min())].End.max()
    strand = genes[(genes.Gene==temp.Gene.min())].Strand.max()
    
    no_vars.append((nv.split("-t")[0],nv,el,refs,alts,nonsyn,nvars,
                    utr5,utr3,exon,intron,chrom,descr,
                    gstart,gend,strand))
    
no_vars = pd.DataFrame(no_vars)
no_vars.columns = ['Gene','Parent','Expected','Refstop','Altstop','Nonsyn','Nvars',
                  'Utr5','Utr3','Exon','Intron','Chrom','Description',
                  'Start','End','Strand']


# In[6]:


samplespath = sorted(glob.glob('../../GENOTYPE/GENES/*/*.csv.gz'))


# In[7]:


todf = []
for s in samplespath:

    sample = s.split('GENES/')[-1].split('/')[0]
    genep = s.split('/')[-1].split('.csv')[0]
    gene_name = s.split('/')[-1].split('-t26')[0]

    temp = pd.read_csv(s).reset_index(drop=True)
    gene = temp.Gene.min()
    assert sample == temp.Sample.min()
    assert gene == genep
        
    cds = temp[(temp.Type==0)]
        
    ref,alt = makeorf(cds)
        
    el = cds.Pos.unique().shape[0]/3-1
        
    ra = ref.translate(to_stop=True)
    aa = alt.translate(to_stop=True)
        
    rl = len(ra)
    al = len(aa)
        
    sr = ref.translate().count('*')
    sa = alt.translate().count('*')
        
    ns = sa - 1
    for i in range(np.min([len(ra),len(aa)])):
        if ra[i]!=aa[i]:
            ns = ns + 1
            
    nvars = temp[(temp.Ref!=temp.Alt)].shape[0]#temp[(temp.Isvar==1)].shape[0]
    utr3 = temp[(temp.Type==3) & (temp.Ref!=temp.Alt)].shape[0]
    utr5 = temp[(temp.Type==5) & (temp.Ref!=temp.Alt)].shape[0]
    inexon = temp[(temp.Type==0) & (temp.Ref!=temp.Alt)].shape[0]
    inintron = temp[(temp.Type==-1) & (temp.Ref!=temp.Alt)].shape[0]
    
    todf.append((sample,gene_name,gene,
                 el,rl,al,sr,sa,ns,
                 nvars,utr5,utr3,inexon,inintron
                ))


# In[8]:


resdf = pd.DataFrame(todf,
        columns=['Strain','Gene','Parent','Expected',
                 'Ref','Alt','Refstop','Altstop','Nonsyn',
                 'Nvars','Utr5','Utr3','Exon','Intron'
                ])

resdf = resdf.merge(genes[['Gene','Chrom','Description',
                           'Start','End','Strand']])

resdf = pd.concat([resdf,no_vars],axis=0,sort=False)


# In[9]:


chroms = '.'.join([str(c) for c in sorted(resdf.Chrom.unique())])
save_path = '../../NOTES/Bt22xFtc555-1_gene_analysis.%s.csv.gz'%chroms
resdf.to_csv(save_path,index=False)
