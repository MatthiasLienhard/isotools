import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_hex
from scipy.stats import beta,nbinom
import seaborn as sns
import pandas as pd
import numpy as np
import logging
logger=logging.getLogger('isotools')

    
def plot_diff_results(result_table, min_support=3, min_diff=.1, grid_shape=(5,5), splice_types=None):
    '''Plots differential splicing results.

    For the first (e.g. most significant) differential splicing events from result_table
    that pass the checks defined by the parameters, 
    the PSI value of the alternative splicing event, as well as the fitted beta model for the groups, is depicted. 

    :param min_support: Minimum number of samples per group supporting the differential event.
    :param min_diff: Minimum PSI group difference.
    :param grid_shape: Number of rows and columns for the figure. 
    :param splice_type: Only events from the splecified splice_type(s) are depicted. If omitted, all types are selected. 
    :return: figure, axes and list of plotted events'''

    plotted=pd.DataFrame(columns=result_table.columns)
    if isinstance(splice_types, str):
        splice_types=[splice_types]
    f,axs=plt.subplots(*grid_shape)
    axs=axs.flatten()
    x=[i/100 for i in range(101)]
    group_names=[c[:-4] for c in result_table.columns if c.endswith('_PSI')][:2]
    groups={gn:[c[:-10] for c in  result_table.columns if c.endswith(gn+'_total_cov')] for gn in group_names}
    logger.debug('groups: %s',str(groups))
    for idx,row in result_table.iterrows():
        logger.debug(f'plotting {idx}: {row.gene}')
        if splice_types is not None and row.splice_type not in splice_types:
            continue
        if row.gene in set(plotted.gene):
            continue
        params_alt={gn:(row[f'{gn}_PSI'],row[f'{gn}_disp']) for gn in group_names}
        psi_gr={gn:[row[f'{sa}_in_cov']/row[f'{sa}_total_cov'] for sa in gr if row[f'{sa}_total_cov']>0] for gn,gr in groups.items()}
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
            max_i=int(params_alt[gn][0]*(len(x)-1))
            ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
            if params_alt[gn][1]>0:
                m,v=params_alt[gn]
                params=((-m*(m**2-m+v))/v,((m-1)*(m**2-m+v))/v)  
                y=beta(*params).pdf(x)
                y[max_i]=beta(*params).pdf(params_alt[gn][0])                
            else:
                y=np.zeros(len(x))
                y[max_i]=1#point mass
            ax2.plot(x,y , color=f'C{i}')
            ax2.tick_params( right=False, labelright=False)
        ax.set_title(f'{row.gene} {row.splice_type}\nFDR={row.padj:.5f}')
        plotted=plotted.append(row)
        if len(plotted)==len(axs):
            break
    f.tight_layout()
    return f,axs,plotted

def plot_embedding(splice_bubbles, method='PCA',prior_count=3, 
        top_var=500,min_total=100,min_alt_fraction=.1,plot_components=[1,2], 
        splice_types='all', labels=True,groups=None,colors=None, ax=None,**kwargs):
    ''' Plots embedding of alternative splicing events. 

    Alternative splicing events are soreted by variance and only the top variable events are used for the embedding. 
    A prior weight is added to all samples proportional to the average fraction of the alternatives, 
    in order to bias poorly covered samples towards the mean and limit their potential to disturb the analysis.

    :param splice_bubbles: The splice bubble table, produced by Transcriptome.alternative_splicing_events().
    :param method: The embedding method, either "PCA" or "UMAP". 
    :param prior_count: Number of prior reads which are added to each sample proportional to the average fraction of the alternatives. 
    :param top_var: Number of alternative splicing events which are used for the embedding. 
    :param min_total: Minimum total coverage over all selected samples.
    :param min_alt_fraction: Minimum fraction of reads supporting the alternative (for both groups combined). 
    :param plot_components: The dimentions to plot (E.g. the components of the PCA)
    :param splice_types: Restrict the analysis on specified splicing event(s).
    :param labels: If True, sample names are printed in the plot next to the corresponding points. 
    :param groups: Set a group definition (e.g. by isoseq.Transcirptome.groups()) to color the datapoints. 
        All samples within one group are depicted. 
    :param colors: List or dict of colors for the groups, if ommited colors are generated automatically.
    :param ax: The axis for plotting.
    :param \**kwargs: Additional keyword parameters are passed to PCA() or UMAP().
    :return: A dataframe with the proportions of the alternative events, the transformed data and the embedding object.'''


    assert method in ['PCA', 'UMAP'], 'method must be PCA or UMAP'
    if method=='UMAP':
        from umap import UMAP as Embedding# pylint: disable-msg=E0611
    else:
        from sklearn.decomposition import PCA as Embedding


    plot_components=np.array(plot_components)
    if isinstance(splice_types, str):
        splice_types=[splice_types]
    if 'all' not in splice_types:        
        splice_bubbles=splice_bubbles.loc[splice_bubbles['splice_type'].isin( splice_types)]
    k=splice_bubbles[[c for c in splice_bubbles.columns if c.endswith('_in_cov')]]
    n=splice_bubbles[[c for c in splice_bubbles.columns if c.endswith('_total_cov')]]
    n.columns=[c[:-10] for c in n.columns]
    k.columns=[c[:-7] for c in k.columns]
    samples=list(n.columns)
    assert all(c1==c2 for c1,c2 in zip (n.columns, k.columns)), 'issue with sample naming of splice bubble table'

    # select samples and assing colors
    if groups is None:
        groups={'all samples':samples}
    else:
        sa_group={sa:gn for gn,salist in groups.items() for sa in salist if sa in samples}
        if len(samples)>len(sa_group):
            samples=[sa for sa in samples if sa in sa_group]
            logger.info('restricting embedding on samples '+', '.join(samples))
            n=n[samples]
            k=k[samples]
    if colors is None:
        cm = plt.get_cmap('gist_rainbow')
        colors={gn:to_hex(cm(i/len(groups))) for i,gn in enumerate(groups)}
    elif isinstance(colors,dict):
        assert all( gn in colors for gn in groups), 'not all groups have colors'
        assert all(is_color_like(c) for c in colors.values()), 'invalid colors'
    elif len(colors)>= len(groups):
        assert all(is_color_like(c) for c in colors), 'invalid colors'
        colors={gn:colors[i] for i,gn in enumerate(groups)}
    else:
        raise ValueError(f'number of colors ({len(colors)}) does not match number of groups ({len(groups)})')
    nsum=n.sum(1)
    ksum=k.sum(1)
    covered=(nsum>=min_total) & (min_alt_fraction < ksum/nsum) & (ksum/nsum < 1-min_alt_fraction)
    n=n.loc[covered]
    k=k.loc[covered]
    #compute the proportions
    scaled_mean=k.sum(1)/n.sum(1)*prior_count
    p=((k.values+scaled_mean[:,np.newaxis])/(n.values+prior_count)).T
    topvar=p[:,p.var(0).argsort()[-top_var:]] # sort from low to high var

    #compute embedding
    kwargs.setdefault('n_components',max(plot_components))
    assert kwargs['n_components']>=max(plot_components), 'n_components is smaller than the largest selected component'
    
        # Linear dimensionality reduction using Singular Value Decomposition of the data to project it to a lower dimensional space. 
        # The input data is centered but not scaled for each feature before applying the SVD.
    embedding=Embedding(**kwargs).fit(topvar)
    axparams=dict(title= f'{method} ({",".join(splice_types)})')
    if method=='PCA':
            axparams['xlabel'] =f'PC{plot_components[0]} ({embedding.explained_variance_ratio_[plot_components[0]-1]*100:.2f} %)' 
            axparams['ylabel'] =f'PC{plot_components[1]} ({embedding.explained_variance_ratio_[plot_components[1]-1]*100:.2f} %)' 
    transformed=pd.DataFrame(embedding.transform(topvar), index=samples)


    if ax is None:
        _,ax=plt.subplots()
    for gr,sa in groups.items():
        ax.scatter(
            transformed.loc[sa,plot_components[0]-1],
            transformed.loc[sa,plot_components[1]-1],
            c=colors[gr], label=gr)
    ax.set(**axparams)
    if labels:
        for idx,(x,y) in transformed[plot_components-1].iterrows():
            ax.text(x,y,s=idx)
    return pd.DataFrame(p.T,columns=samples, index=k.index ),transformed, embedding
  
#plots
def plot_bar(df,ax=None, drop_categories=None, legend=True, annotate=True,rot=90, bar_width=.5,**axparams):   
    '''Depicts data as a barplot.

    This function is intended to be called with the result from 
    isoseq.Transcriptome.filter_stats() or isoseq.Transcriptome.altsplice_stats().
    
    :param df: Pandas dataframe with the data to plot. 
    :param ax: the axis for the plot.
    :param drop_categories: Specify columns from df to drop. 
    :param legend: If True, add a legend.
    :param annotate: If True, print numbers / fractions in the bars.
    :param rot: Set rotation of the lables.
    :param bar_width: Set relative width of the plotted bars.
    :param \**axparams: Additional keyword parameters are passed to ax.set().'''

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
    fractions.drop(dcat).plot.bar( ax=ax, legend=legend, width=bar_width, rot=rot)   
    #add numbers 
    if annotate:
        numbers=[int(v) for c in df.drop(dcat).T.values for v in c]
        frac=[v for c in fractions.drop(dcat).T.values for v in c]
        for n,f,p in zip(numbers,frac,ax.patches):
            small=f<max(frac)/2
            #contrast=tuple(1-cv for cv in p.get_facecolor()[:3])
            contrast='white' if np.mean(p.get_facecolor()[:3])<.5 else 'black'
            ax.annotate(f' {f/100:.2%} ({n}) ', (p.get_x()+p.get_width()/2 , p.get_height() ) ,ha='center',  va='bottom' if small else 'top' ,rotation=90, color='black' if small else contrast, fontweight='bold')
    ax.set(**axparams)
    
    return ax

def plot_distr(counts,ax=None,density=False,smooth=None,  legend=True,fill=True,**axparams):
    '''Depicts data as density plot.

    This function is intended to be called with the result from 
    isoseq.Transcriptome.transcript_length_hist(), isoseq.Transcriptome.transcripts_per_gene_hist(), 
    isoseq.Transcriptome.exons_per_transcript_hist(), isoseq.Transcriptome.downstream_a_hist(), 
    isoseq.Transcriptome.direct_repeat_hist() or isoseq.Transcriptome.transcript_coverage_hist().
    
    :param df: Pandas dataframe with the data to plot. 
    :param ax: The axis for the plot.
    :param density: Scale the data by the total. 
    :param smooth: Ews smoothing span.
    :param legend: If True, add a legend.
    :param fill: If set, the area below the lines are filled with half transparent color. 
    :param \**axparams: Additional keyword parameters are passed to ax.set().'''
    #maybe add smoothing
    x=[sum(bin)/2 for bin in counts.index]
    sz=[bin[1]-bin[0] for bin in counts.index]
    if ax is None:
        _, ax=  plt.subplots()
    if density: 
        counts=(counts/counts.sum())
        if 'ylabel' in axparams and 'density' not in axparams['ylabel']:
            axparams['ylabel']+=' density'
        else:
            axparams['ylabel']='density'
    else:
        axparams.setdefault('ylabel','# transcripts')
    if smooth:
        counts=counts.ewm(span=smooth).mean()
    for gn,gc in counts.items():
        ax.plot(x, gc/sz, label=gn)    
        if fill:
            ax.fill_between(x, 0, gc/sz,  alpha=.5)
    #ax.plot(x, counts.divide(sz, axis=0))
    ax.set(**axparams)
    if legend:
        ax.legend()
    return ax

def plot_saturation(isoseq=None,ax=None,cov_th=2,expr_th=[.5,1,2,5,10],x_range=(1e4,1e7,1e4),legend=True,label=True, **axparams):
    '''Plots Negative Binomial model to analyze the saturation of LRTS data.
    
    Saturation (e.g. the probability to observe a transcript of interest in the sample) is dependent on the sequencing depth (number of reads), 
    the concentration of the transcripts of interest in the sample (in TPM),
    and the requested coverage of the transcript in the data (minimum number of reads per transcript). 
    This function models the relation with a Negative Binomial distribution, to help estimate the required sequencing depth.

    :param isoseq: If provided, the sequencing depth of samples from this isotools.Transcriptome object are depicted as vertical lines. 
    :param ax: The axis for the plot. 
    :param cov_th: The requested coverage, e.g. the minimum number of reads per transcript.
    :param expr_th: A list of transcript concentrations in TPM for transcripts of interest. 
    :param x_range: Specify the range of the x axis (e.g. the sequencing depth)
    :param legend: If set True, a legend is added to the plot. 
    :param label: If set True, the sample names and sequencing depth from the isoseq parameter is printed in the plot. 
    :param \**axparams: Additional keyword parameters are passed to ax.set().'''
    if ax is None:
        _, ax=  plt.subplots()
    k=np.arange(*x_range)
    axparams.setdefault('title','Saturation Analysis') #[nr],{'fontsize':20}, loc='left', pad=10)
    axparams.setdefault('ylabel',(f'probaility of sampling at least {cov_th} transcript{"s" if cov_th>1 else ""}'))
    axparams.setdefault('ylim',(0,1))
    axparams.setdefault('xlabel','number of reads [million]')
    n_reads=isoseq.sample_table.set_index('name')['nonchimeric_reads'] if isoseq is not None else {}
    for tpm_th in expr_th:
        chance = nbinom.cdf(k-cov_th, n=cov_th, p=tpm_th*1e-6) # 0 to k-cov_th failiors
        ax.plot(k/1e6, chance, label=f'{tpm_th} TPM')
    for sa,cov in n_reads.items():
        ax.axvline(cov/1e6, color='grey', linestyle='--')
        if label:
            ax.text((cov+(k[-1]-k[0])/200)/1e6,0.1,f'{sa} ({cov/1e6:.2f} M)',rotation=-90)
    ax.set(**axparams)
    
    if legend:
        ax.legend()
    return ax