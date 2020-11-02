import matplotlib.colors as plt_col
import matplotlib.pyplot as plt
from scipy.stats import beta

import seaborn as sns
import pandas as pd
import numpy as np
from math import log10
import logging
logger=logging.getLogger('isotools')

    
def plot_diff_results(result_table, min_support=3, min_diff=.1, grid_shape=(5,5), splice_types=None):
    plotted=pd.DataFrame(columns=result_table.columns)
    if isinstance(splice_types, str):
        splice_types=[splice_types]
    f,axs=plt.subplots(*grid_shape)
    axs=axs.flatten()
    x=[i/1000 for i in range(1001)]
    group_names=[c[:-9] for c in result_table.columns if c[-9:]=='_fraction'][:2]
    groups={gn:[c[4:] for c in  result_table.columns if c[:4]=='cov_' and c.endswith(gn)] for gn in group_names}
    logger.debug('groups: %s',str(groups))
    for idx,row in result_table.iterrows():
        logger.debug(f'plotting {idx}: {row.gene}')
        if splice_types is not None and row.splice_type not in splice_types:
            continue
        if row.gene in set(plotted.gene):
            continue
        params_alt={gn:(row[f'{gn}_fraction'],row[f'{gn}_disp']) for gn in group_names}
        psi_gr={gn:[row[f'cov_{sa}']/row[f'span_cov_{sa}'] for sa in gr if row[f'span_cov_{sa}']>0] for gn,gr in groups.items()}
        support={s: sum(abs(i-params_alt[s][0]) < abs(i-params_alt[o][0]) for i in psi_gr[s]) for s,o in zip(group_names, reversed(group_names))}
        if any(sup<min_support for sup in support.values()):
            logger.debug(f'skipping {row.gene} with {support} supporters')
            continue
        if abs(params_alt[group_names[0]][0]-params_alt[group_names[1]][0])<min_diff:  
            logger.debug(f'{row.gene} with {"vs".join(str(p[0]) for p in params_alt.values())}')  
            continue
        #get the paramters for the beta distiribution
        #print(param)
        ax=axs[len(plotted)]
        #ax.boxplot([mut,wt], labels=['mut','wt'])
        sns.swarmplot(data=pd.DataFrame(list(psi_gr.values()), index=psi_gr).T, ax=ax, orient='h')
        for i,gn in enumerate(group_names):
            if params_alt[gn][1]>0:
                m,v=params_alt[gn]
                params=((-m*(m**2-m+v))/v,((m-1)*(m**2-m+v))/v)  
                y=beta(*params).pdf(x)
                ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
            else:
                y=np.zeros(len(x))
                y[int(params_alt[gn][0]*(len(x)-1))]=1
            ax2.plot(x,y , color=f'C{i}')
            ax2.tick_params( right=False, labelright=False)
        ax.set_title(f'{row.gene} {row.splice_type}\nFDR={row.padj:.5f}')
        plotted=plotted.append(row)
        if len(plotted)==len(axs):
            break
    f.tight_layout()
    return f,axs,plotted

def plot_embedding(embedding, col=None, groups=None, labels=True, ax=None, comp=[0,1]):    
    if groups is not None:
        if col is None:
            cm = plt.get_cmap('gist_rainbow')
            col=[cm(i/len(groups)) for i in range(len(groups))]
    if ax is None:
        _,ax=plt.subplots()
    ax.scatter(
        embedding[comp[0]],
        embedding[comp[1]],
        c=col)
    if labels:
        for idx,(x,y) in embedding[comp].iterrows():
            ax.text(x,y,s=idx)
  
#plots
def plot_bar(df,ax=None, drop_categories=None, legend=True,**axparams):    
    if ax is None:
        _, ax=  plt.subplots()
    if 'total' in df.index:
        total=df.loc['total']
        df=df.drop('total')
    else:
        total=df.sum()
    fractions=(df/total*100)
    if drop_categories is None:
        dcat=[]
    else:
        dcat=[d for d in drop_categories if d in df.index]
    fractions.drop(dcat).plot.bar( ax=ax)   
    #add numbers 
    numbers=[int(v) for c in df.drop(dcat).T.values for v in c]
    frac=[v for c in fractions.drop(dcat).T.values for v in c]
    for n,f,p in zip(numbers,frac,ax.patches):
        small=f<max(frac)/2
        #contrast=tuple(1-cv for cv in p.get_facecolor()[:3])
        contrast='white' if np.mean(p.get_facecolor()[:3])<.5 else 'black'
        ax.annotate(f' {f/100:.2%} ({n}) ', (p.get_x()+p.get_width()/2 , p.get_height() ) ,ha='center',  va='bottom' if small else 'top' ,rotation=90, color='black' if small else contrast, fontweight='bold')
    ax.set(**axparams)
    if legend:
        ax.legend()
    return ax

def plot_distr(counts,ax=None,density=False,smooth=None,  legend=True,fill=True,**axparams):
    #maybe add smoothing
    x=[sum(bin)/2 for bin in counts.index]
    sz=[bin[1]-bin[0] for bin in counts.index]
    if ax is None:
        _, ax=  plt.subplots()
    if density: 
        counts=(counts/counts.sum())
    if smooth:
        counts=counts.ewm(span=smooth).mean()
    if fill:
        for gn,gc in counts.items():
            ax.fill_between(x, 0, gc/sz, label=gn, alpha=.5)
    ax.plot(x, counts.divide(sz, axis=0))
    ax.set(**axparams)
    if legend:
        ax.legend()
    return ax
