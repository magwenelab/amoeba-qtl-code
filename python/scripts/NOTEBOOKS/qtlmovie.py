import numpy as np, pandas as pd, scipy.stats as ss, qtlfunctions as qtlfun, os
from matplotlib import pyplot as plt


def plot_pheno(ax,plotbox,gxp,
               xjitter=None,
               yjitter = None,
               plotreg=True,
               ylims = (0,60),
               xlims = (-0.5,1.5),
               swarmcolor='k',
               swarmalpha=0.1,
               regcolor='r',
               xticklabels = [],
               yticklabels = [],
               yticks = [],
               ylabel = None,
               colors = ['tab:blue','tab:orange'],
               annot = None,
               fs = 10):
    """
    Plot phenotype by genotype swarm
    """
    if yjitter is None:
        yjitter = np.zeros(gxp.shape[0])
    if xjitter is None:
        xjitter = np.zeros(gxp.shape[0])
    if len(yticks) > 0 and len(yticklabels) == 0:
      yticklabels = yticks
        
    pax = ax.inset_axes(plotbox, transform=ax.transData)
    pax.set_xlim(*xlims)
    pax.set_ylim(*ylims)
    pax.set_xticklabels(xticklabels,fontsize=fs)
    pax.set_yticks(yticks)
    pax.set_yticklabels(yticklabels[:len(yticks)],fontsize=fs)
    pax.set_ylabel(ylabel,fontsize=fs-4)
    
    pax.plot(gxp.Genotype.values+xjitter,
             gxp.Phenotype.values+yjitter,
             '.',rasterized=True,
             color=swarmcolor,alpha=swarmalpha)
    
    if plotreg:
        pax.plot(gxp.groupby('Genotype').mean().index,
             gxp.groupby('Genotype').mean().Phenotype.values,
             color=regcolor,alpha=0.7,linewidth=1)
        
        for g in [0,1]:
            pax.plot(g, gxp[(gxp.Genotype==g)].Phenotype.mean(),
                     '^',color=colors[int(g)],alpha=0.85)
    
    for g in [0,1]:
        yeast = pax.inset_axes([0.15 if g ==0 else 0.65,-0.2,.3,.2],
                               transform=pax.transAxes)
        yeast.plot(0,0,'o',color=colors[g],ms=7)
        yeast.plot(0.8,0.4,'o',color=colors[g],ms=3)
        yeast.set_xlim(-1,1.5)
        yeast.set_ylim(-1,1.5)
        yeast.axis('off')
        yeast.patch.set_alpha(0)
        
    if annot is not None:
        pax.set_title(annot,fontsize=fs)
    pax.patch.set_alpha(0)
    return pax


def QTLbin(df,qtln,grouper='Chrom',ns=1000,qtlns=200):
    """
    Bins a qtl dataframe (DF) for iterating QTL figures.
    
    Chromosomes with a QTL present should be listed in QTLN.
    The dataframe (DF) is grouped by the variable in GROUPER; default is 'Chrom'
    NS and QTLNS are number of row indexs within each chromosome bin 
    across the genome or for chromosomes with QTLs (listed in QTLN).
    
    Returns a list of dataframe indices.
    """
    kix = []
    for i,c in df.groupby(grouper):
        ixn = ns if i not in qtln else qtlns
        kc = c.index[::ixn]
        kix.append(kc)
    return list(np.concatenate(kix)) + [df.index.max()]

def clearfigs(path,figtype='png'):
    """
    Removes figures (ending in FIGTYPE (ie .png)) along a specified PATH
    """
    newpath = '%s*%s'%(path,figtype)
    
    w = os.system('rm %s'%newpath)
    if w != 0:
        print('Failed to remove figures on specified path:\n%s'%path)

        
def recombprog(geno,mappop,nrecsegs):
    """
    Given a genotype dataframe (GENO) and a list of segregants (MAPPOP)
    returns a set number (NRECSEGS) of segregant names nearing the 
    average number crossovers.
    """
    recsegs = geno[mappop].diff(axis=0).abs().sum().sort_values().index.values[
        int((len(mappop)- nrecsegs)/2):int((len(mappop) + nrecsegs)/2)]
    return recsegs[:nrecsegs]
    
        
def QTLmover(mpvaldf,Mappop,pheno,QTLixs,QTLn,relative_path,
              max_pval=0,yaxis_max = None,
              recsegs=None,chrommap=None,kix=None,
              phenomax=None,phenomod=20,phenomin=0,
              x_jitter=None,y_jitter=None,x_std = 0.1,y_std=0.1,
              pval_ymod=13,xwidth=11**6,qlast=-1,
              pval_tick_mod=20,yheight_mod=0.75,
              grouper='Chrom',ns=1000,qtlns=200,
              kcmap = ['k','lightgrey','tab:grey'],
              myfs=10,
              mydpi=100,
              ylabely=0.35,
              qtlname='QTL',
              annot=True,
              pheno_ylabel = None,
              haplo_xlabel='Chromosome',
              haplo_ylabel='Example F$_1$\nHaplotypes\n',
              colors = ['tab:blue','tab:orange'],
              poscol='Cumpos',
              R2 = 'R$^2$ = %s',
              savetype='png',
              facecolor='w',
              figsize=(12,5),
              dnamod=0.075,
              chrompad = 200000,
              xpad = 700000,
              nsegs = 8,
              qtl = 'QTL',
              pheno_rounder = -1,
              pval_rounder = -1,
              yheight_rounder = -1,
             ):
    """
    Given a QTL dataframe (MPVALDF) iteratively plot manhattan plots, haplotype maps, 
    and genotype by phenotype plots across the genome for use in QTL GIFs.
    """
    
    togif = []
    qtlgif = []
    QTLplots = []
    
    if phenomax is None:
        phenomax = np.round(np.max(pheno),pheno_rounder)

    if max_pval == 0:
        max_pval = np.round(mpvaldf[qtl].max(),pval_rounder)

    xaxismin = -2*xpad
    yheight = np.round(yheight_mod*max_pval,yheight_rounder)
    pvalu_ticks = np.arange(0,np.round(max_pval)+pval_tick_mod,pval_tick_mod)
    phenoticks = np.arange(phenomin,phenomax+phenomod,phenomod)

    ## Set jitter of x and y axis
    if x_jitter is None:
        x_jitter = np.random.normal(0,x_std,len(pheno))
    
    if y_jitter is None:
        y_jitter = np.random.normal(0,y_std,len(pheno))
        
    if kix is None:
        kix = QTLbin(mpvaldf,QTLn,grouper=grouper,ns=ns,qtlns=qtlns)
        
    if chrommap is None:
        chrommap = qtlfun.chrommap(mpvaldf)
        chrompad = np.cumsum(chrompad*np.ones(len(chrommap)))
        chrommap['Cumlen'] = chrommap.Cumlen + chrompad
        chrommap['Midpts'] = chrommap.Midpts + chrompad

    if yaxis_max is None:
      yaxis_max = 2*max_pval+pval_ymod
        
    maxx = chrommap.Cumlen.max()+chrommap.Length.max()
        
    if poscol not in mpvaldf.columns.tolist():
        mpvaldf = mpvaldf.merge(chrommap)
        mpvaldf[poscol] = mpvaldf.Pos + mpvaldf.Cumlen
        
    qtlixs = mpvaldf.index.values

    if recsegs is None:
        recsegs = recombprog(mpvaldf,Mappop,nsegs)
        
    for j,i in enumerate(kix+[np.max(kix)]):
    
    ## Make a figure and its subplots via grid spec
        fig3 = plt.figure(figsize=figsize)
        fig3.set_facecolor(facecolor)
        gs = fig3.add_gridspec(3, 3)
        f1_ax1 = fig3.add_subplot(gs[:2, :])
        f2_ax1 = fig3.add_subplot(gs[2,:])
    
    ## Set save path
        savepath=relative_path+'%s.%s'%(str(j).zfill(len(str(len(kix)))+1),savetype)
        togif.append(savepath)

    ## Gather index of QTL values to plot
        tocolor = qtlixs[(qtlixs<=i)]
    
    ## Gather the current and last index of QTL values
        thisix = tocolor[(tocolor>=qlast)]
        qlast = np.max(tocolor) 
    
    ## Gather positional info
        chrom = mpvaldf.loc[thisix]
        posmin = mpvaldf.loc[thisix,poscol].min()
        posmax = mpvaldf.loc[thisix,poscol].max()
        posmean = np.mean([posmin,posmax])
        pvalmax = mpvaldf.loc[thisix,qtlname].max()
    
    ## Plot the QTL
        qtlfun.manhattan(mpvaldf.loc[tocolor,:],chrmap=chrommap,ax = f1_ax1,
                          ylabel=None,xlabel=None,fs = myfs,mydpi=mydpi,cmap=kcmap)
    
        plt.sca(f1_ax1)
        plt.xlim(xaxismin,maxx)
        plt.xticks(chrommap.Midpts,[])
        plt.ylim(-0.1*max_pval,yaxis_max)
        plt.yticks(pvalu_ticks)
        plt.vlines(xaxismin,0,np.max(pvalu_ticks),color='k')
        plt.ylabel(qtlfun.pvalylable,y=ylabely)
    
    ## Plot phenotype by genotype
        plotbox = [int(posmax-(xwidth/2)),
               pval_ymod+int(pvalmax),
               xwidth,yheight]
    
        thisix_max = mpvaldf.loc[thisix,qtlname].idxmax()
        this_gxp = qtlfun.qtldf(mpvaldf,pheno,Mappop,qtlix=thisix_max)
    
        if j < len(kix):
            pax = plot_pheno(f1_ax1,plotbox,this_gxp,x_jitter,
              ylims=(phenomin,phenomax),ylabel=None)#=pheno_ylabel if thisix_max in QTLixs else None)
            [pax.spines[ts].set_visible(False)for ts in ['left','top','right']]
            pax.set_yticks([])
            plt.vlines(posmax,pvalmax,pval_ymod+int(pvalmax),color='k',linewidth=0.8)
    
        else :
            pass
    
        if thisix_max in QTLixs:
        
            [pax.spines[ts].set_visible(True) for ts in ['left']]
            pax.set_yticks(phenoticks)

            QTLplots.append((plotbox,this_gxp))
            qtlgif.append(savepath)
        
        for (qtlbox,qtlgxp) in QTLplots:
        
            if annot:
                r2 = R2%str(np.round(ss.pearsonr(qtlgxp.Genotype.values,
                                  qtlgxp.Phenotype.values)[0]**2,2))
            else:
                r2 = None
        
            qtlax = plot_pheno(f1_ax1,qtlbox,qtlgxp,x_jitter,colors=colors,
                           ylims=(phenomin,phenomax),annot=r2,yticks=phenoticks,
                           ylabel=pheno_ylabel)
            qtlax.set_yticks(phenoticks)
            [qtlax.spines[ts].set_visible(False) for ts in ['top','right']]

    ## Plot haplotypes
        plt.sca(f2_ax1)
        plt.xlim(xaxismin,maxx)
        plt.xticks(chrommap.Midpts,chrommap.index.values+1,fontsize=myfs)
        plt.yticks([])
        plt.xlabel(haplo_xlabel,fontsize=myfs)
        plt.ylabel(haplo_ylabel,fontsize=myfs)

        for si,s in enumerate(recsegs):
            stemp = mpvaldf[[s,poscol]]

            for gi in [0,1]:
                symod = (dnamod) if gi == 0 else (-dnamod),
                plt.plot(stemp[(stemp[s]==gi)].Cumpos.values,
                         si*np.ones(len(stemp[(stemp[s]==gi)][s].values)
                        ) + symod*np.ones(len(stemp[(stemp[s]==gi)][s].values)),
                         '.',ms=1, alpha=0.5, rasterized=True,
                         color= colors[int(gi)])
            
            sgt = this_gxp.loc[s].Genotype    
            sgtcolor = colors[int(sgt)]
            
            ## Plot the cartoon yeast
            plt.plot(0-xpad/2,si,'o',color=sgtcolor,ms=6,rasterized=True)
            plt.plot(0-xpad/2+(xpad/5),si+0.1,'o',color=sgtcolor,ms=3,rasterized=True)
    
    ## Turn of spins and adjust plots
        [f1_ax1.spines[ts].set_visible(False) for ts in ['top','left','right','bottom']]
        [f2_ax1.spines[ts].set_visible(False) for ts in ['top','left','right']]
        plt.subplots_adjust(hspace=0.0)  
        plt.savefig(savepath,dpi=mydpi,bbox_inches='tight')
        plt.close()
        
    return togif,qtlgif,recsegs


def QTLgif(path,images,qtlimages):
    """
    Formats and submits an image majics command converting images to a GIF.
    """
    pass