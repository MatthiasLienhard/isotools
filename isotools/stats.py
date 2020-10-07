from math import log10
import matplotlib.colors as plt_col
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import pandas as pd
import numpy as np
from pysam import  AlignmentFile
from scipy.stats import binom,norm, chi2, betabinom, beta
from scipy.special import gammaln
from scipy.optimize import minimize
import statsmodels.stats.multitest as multi
import isotools.splice_graph
from isotools.transcriptome import junctions_from_cigar, Gene
from tqdm import tqdm
import logging
from heapq import nlargest
from umap import UMAP
from sklearn.decomposition import PCA
import itertools

log=logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s: %(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.handlers=[]
log.addHandler(log_stream)

def overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False

def embedding(isoseq, genes=None, reducer='PCA',n_top_var=500, n_subsample=None,samples=None, min_cov=20):
    ''' Compute embedding of most variable junction for each gene if genes are provided as a list (e.g. genes=[gene_name1, gene_name2 ...]). 
        Alternatively, if genes is a dict, junctions can be provided as a list (e.g. genes[gene_name]=[(jstart,jend)]).
        If genes is None, all genes with coverage > min_cov are considered. 
    '''
    joi=dict()
    if samples is not None:
        s_filter=[True if sn in samples else False for sn in isoseq.runs]
    else:
        s_filter=[True]*len(isoseq.runs)
    if genes is None:        
        genes=[g for i,g in enumerate(isoseq) if all(g.splice_graph.weights[s_filter].sum(1)>=min_cov)] #all covered genes
        log.info(f'found {len(genes)} genes > {min_cov} reads in all {sum(s_filter)} selected samples.')
    else:
        n_in=len(genes)
        if isinstance(genes, dict):
            joi=genes
        genes=[isoseq[g] for i,g in enumerate(genes) if g in isoseq and all(isoseq[g].splice_graph.weights[s_filter].sum(1)>=min_cov)] #all covered genes
        log.info(f'{len(genes)} of {n_in} genes have > {min_cov} reads in all {sum(s_filter)} selected samples.')
    if n_subsample is not None and n_subsample<len(genes):
        genes=[genes[i] for i in np.random.choice(range(len(genes)), n_subsample)]
        log.info(f'sampled {n_subsample} random genes')        
    data={}
    # since junctions of a single gene are often correlated, I select one top variable junction per gene
    for g in tqdm(genes):   
        sg=g.splice_graph
        gene_coverage=sg.weights[s_filter].sum(1)
        if any(cov<min_cov for cov in gene_coverage):
            continue
        if g.name in joi and joi[g.name]:#
            for j in joi[g.name]:
                start=next(n for n in sg if n.end==j[0])            
                jcov=sum(sg.weights[s_filter,trid] for trid,suc in start.suc.items() if sg[suc].start==j[1])/gene_coverage
                data[f'{g.name}_{j[0]}_{j[1]}']=jcov
        else:
            jcov={}                 
            for n in sg:
                #print(n)
                for trid,suc in n[3].items():
                    #if all(sg.weights[s_filter,trid] < min_cov):
                    #    continue
                    j_name=f'{g.name}_{n[1]}_{sg[suc][0]}'
                    jcov.setdefault(j_name,np.zeros(len(gene_coverage)))
                    jcov[j_name]+=(sg.weights[s_filter,trid]) /gene_coverage
            if jcov:
                top_idx=pd.DataFrame(jcov).var(0).idxmax()
                data[top_idx]=jcov[top_idx]
    data=pd.DataFrame(data,index=isoseq.infos['sample_table'].loc[s_filter,'name'])
    #filter for highest variance
    log.info(f'selecting the most variable junction from {min(data.shape[1],n_top_var)} genes.')
    topvar=data[data.var(0).nlargest(n_top_var).index]

    #compute embedding
    if reducer=='PCA':
        # Linear dimensionality reduction using Singular Value Decomposition of the data to project it to a lower dimensional space. 
        # The input data is centered but not scaled for each feature before applying the SVD.
        reducer=PCA()
    elif reducer=='UMAP':
        reducer==UMAP()
    elif reducer is None:
        return topvar    
    return pd.DataFrame(reducer.fit_transform(topvar), index=topvar.index)

def plot_embedding(embedding, col=None, groups=None, labels=True, ax=None):    
    if groups is not None:
        if col is None:
            cm = plt.get_cmap('gist_rainbow')
            col=[cm(i/len(groups)) for i in range(len(groups))]
    if ax is None:
        _,ax=plt.subplots()
    ax.scatter(
        embedding.loc[:, 0],
        embedding.loc[:, 1],
        c=col)
    if labels:
        for idx,(x,y) in embedding.loc[:,:1].iterrows():
            ax.text(x,y,s=idx)
    


#sashimi plots
def sashimi_plot_bam(g,ax=None,text_width=.02, text_height=1, title=None,group=None, high_cov_th=.1,chrom_map=None,junctions_of_interest=None, exon_color='blue', 
                low_cov_junctions={'color':'grey','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':1,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':2,'draw_label':True}):
    jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    delta=np.zeros(g.end-g.start)
    chrom=g.chrom
    if title is None:
        title=g.name    
    if text_width<1:
        text_width=(g.end-g.start)*text_width
    #cov=np.array([sum(illu[pos] for i, illu in enumerate(g.illumina_coverage) if groups is None or i in groups) for pos in range(g.start, g.end)])
    cov=np.zeros(g.end-g.start)
    junctions={}
    for i,illu in enumerate(g.illumina_coverage):
        if group is None or i in group:
            cov+=illu.profile
            for k,v in illu.junctions.items():            
                junctions[k]=junctions.get(k,0)+v     
                   
    total_weight=max(cov)
    if high_cov_th<1:
        high_cov_th=high_cov_th*total_weight

    if ax is None:
        fig,ax = plt.subplots(1)
        
    ax.fill_between(range(g.start, g.end), 0, np.log10(cov, where=cov>0, out=np.nan*cov),facecolor=exon_color )
    
    textpositions=[]
    for (x1,x2),w in junctions.items():
        if x1<=g.start or x2>=g.end:
            #todo: this seems to happen at some chrMT genes?
            log.debug(f'attempt to plot junction ({(x1,x2)}) outside of gene: {g.__str__()}')
            continue
        y1=cov[x1-g.start-1]+.5
        y2=cov[x2-g.start]+.5
        center=(x1+x2)/2 
        width=x2-x1
        bow_height=text_height
        while any(overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
            bow_height+=text_height
        textpositions.append((center, log10(max(y1,y2))+bow_height))
        if y1<y2: 
             bow_height=(log10(y2/y1)+bow_height,bow_height)
        elif y1>y2:
             bow_height=(bow_height,bow_height+log10(y1/y2))
        else:
            bow_height=(bow_height,bow_height)
        if junctions_of_interest is not None and (x1,x2) in junctions_of_interest:
            priority=2
        elif w< high_cov_th:
            priority=0
        else:
            priority=1
        bow1=patches.Arc((center, log10(y1)), width=width, height=bow_height[0]*2,theta1=90, theta2=180,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        bow2=patches.Arc((center, log10(y2)), width=width, height=bow_height[1]*2,theta1=0, theta2=90,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if jparams[priority]['draw_label']:
            txt=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,w,horizontalalignment='center', verticalalignment='bottom',bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
        #bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

    
    ax.set_xlim(g.start-100, g.end+100)
    if textpositions:
        ax.set_ylim(-text_height,max(tp[1] for tp in textpositions)+2*text_height)
    else:
        ax.set_ylim(-text_height,3) #todo: adjust y axis and ticklabels to coverage
    ax.set(frame_on=False)
    ax.set_yticks([0,1,2,3])
    ax.set_yticklabels([1,10,100,1000])
    #ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    ax.set_title(title)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    return(ax)

def sashimi_plot(g, ax=None,text_width=.02, arc_type='coverage',text_height=1,
                title=None,group=None,  high_cov_th=.1,junctions_of_interest=None,  exon_color='blue', 
                low_cov_junctions={'color':'grey','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':1,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':2,'draw_label':True}):
    jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    
    sg=g.splice_graph
    if title is None:
        title=g.name
    if group is None:
        group=list(range(sg.weights.shape[0])) #all
    
    boxes=[(node[0], node[1], sg.weights[np.ix_(group,[i for i in set(node[2]).union(node[3])])].sum()) for node in sg]
    if text_width<1:
        text_width=(sg[-1][1]-sg[0][0])*text_width
    total_weight=sg.weights[group,:].sum()
    if high_cov_th<1:
        high_cov_th=high_cov_th*total_weight
    
    #idx=list(range(len(sg)))
    arcs=[]
    for i,(_,ee, _, suc) in enumerate(sg._graph):
        weights={}
        for tr,next_i in suc.items():
            weights.setdefault(next_i,0)
            weights[next_i]+=sg.weights[np.ix_(group,[tr])].sum()
        arcs_new=[(ee,boxes[i][2],sg[next_i][0],boxes[next_i][2],w) for next_i, w in weights.items() if sg[next_i][0]>ee and w>0]
        if arcs_new:
            arcs.extend(arcs_new)
    if ax is None:
        fig,ax = plt.subplots(1)
        
    for st, end, h in boxes:
        if h>0:
            rect = patches.Rectangle((st,0),(end-st),log10(h),linewidth=1,edgecolor=exon_color,facecolor=exon_color,zorder=5)
            ax.add_patch(rect)
    textpositions=[]
    for x1,y1,x2,y2,w in arcs:        
        center=(x1+x2)/2 
        width=x2-x1
        bow_height=text_height
        while any(overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
            bow_height+=text_height
        textpositions.append((center, log10(max(y1,y2))+bow_height))
        if y1<y2: 
            bow_height=(log10(y2/y1)+bow_height,bow_height)
        elif y1>y2:
            bow_height=(bow_height,bow_height+log10(y1/y2))
        else:
            bow_height=(bow_height,bow_height)
        if junctions_of_interest is not None and (x1,x2) in junctions_of_interest:
            priority=2
        elif w< high_cov_th:
            priority=0
        else:
            priority=1
        bow1=patches.Arc((center, log10(y1)), width=width, height=bow_height[0]*2,theta1=90, theta2=180,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        bow2=patches.Arc((center, log10(y2)), width=width, height=bow_height[1]*2,theta1=0, theta2=90,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if arc_type=='coverage':
            lab=str(w)
        else:                #fraction
            lab=f'{w/total_weight:.1%}'
            if arc_type=='both':
                lab=str(w)+' / '+lab
        if jparams[priority]['draw_label']:
            _=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,lab,
                    horizontalalignment='center', verticalalignment='bottom',zorder=10+priority,
                    bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
        #bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

    #ax.set_yscale('log') 
    ax.set_xlim(g.start-100, g.end+100)
    if textpositions:
        ax.set_ylim(-text_height,max(tp[1] for tp in textpositions)+2*text_height)
    else:
        ax.set_ylim(-text_height,3) #todo: adjust y axis and ticklabels to coverage
    ax.set(frame_on=False)
    ax.set_yticks([0,1,2,3])
    ax.set_yticklabels([1,10,100,1000])


    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))

    #ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    #ax.set_xscale(1e-6, 'linear')
    ax.set_title(title)
    return(ax)

def gene_track(genes, ax=None,title=None, draw_exon_numbers=True, color='blue'):
    contrast='white' if np.mean(plt_col.to_rgb(color))<.5 else 'black'
    if ax is None:
        _,ax = plt.subplots(1)    
    if isinstance(genes, Gene):
        genes=[genes]
    start=min(g.start for g in genes)
    end=max(g.end for g in genes)
    assert all(g.chrom==genes[0].chrom for g in genes), "all chromosomes must be equal"
    i=0
    for g in genes:
        for tr in g.transcripts.values():
            #line from TSS to PAS at 0.25
            ax.plot((tr['exons'][0][0], tr['exons'][-1][1]), [i+.25]*2, color=color)
            #idea: draw arrow to mark direction?
            for j,(st, end) in enumerate(tr['exons']):
                if 'CDS' in tr and tr['CDS'][0] <= end and tr['CDS'][1] >= st:#CODING exon
                    c_st,c_end=max(st,tr['CDS'][0]), min(tr['CDS'][1],end) #coding start and coding end
                    if c_st > st: #first noncoding part
                        rect = patches.Rectangle((st,i+.125),(c_st-st),.25,linewidth=1,edgecolor=color,facecolor=color) 
                        ax.add_patch(rect)  
                    if c_end < end: #2nd noncoding part
                        rect = patches.Rectangle((c_end,i+.125),(end-c_end),.25,linewidth=1,edgecolor=color,facecolor=color) 
                        ax.add_patch(rect)  
                    #Coding part                
                    rect = patches.Rectangle((c_st,i),(c_end-c_st),.5,linewidth=1,edgecolor=color,facecolor=color) 
                    ax.add_patch(rect)  
                else: #non coding
                    rect = patches.Rectangle((st,i+.125),(end-st),.25,linewidth=1,edgecolor=color,facecolor=color)
                    ax.add_patch(rect)  
                if draw_exon_numbers:
                    enr=j+1 if g.strand=='+' else len(tr['exons'])-j
                    ax.text((st+end)/2,i+.25,enr,ha='center', va='center', color=contrast).set_clip_on(True)    #bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)
            i+=1
    if title is None:
        title=', '.join(g.name for g in genes)
    ax.set_title(title)
    ax.set(frame_on=False)
    ntr=sum(g.n_transcripts for g in genes)
    ax.set_yticks([i+.25 for i in range(ntr)])
    ax.set_yticklabels(tr for g in genes for tr in g.transcripts)
    ax.tick_params(left=False)
    ax.set_ylim(-.5,ntr+1)
    ax.set_xlim(start-100, end+100)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    return ax


# differential splicing
def proportion_test(x,n):
    # Normal approximation
    #x,n should be lenght 2(the two groups)
    #tests H0: proportions are equal vs H1: proportions are different (two sided)
    x=[xi.sum() for xi in x]
    n=[ni.sum() for ni in n]
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    z=abs(p1[0]-p1[1])/np.sqrt(p0*(1-p0)*(1/n[0]+1/n[1]))
    return(2*norm.sf(z)), (*p1,p0)#two sided alternative
    
def binom_lr_test(x,n):
    # likelihood ratio test
    # x,n should be length 2 (the two groups)
    # principle: log likelihood ratio of M0/M1 is chi2 distributed
    x=[xi.sum() for xi in x]
    n=[ni.sum() for ni in n]
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    # calculate the log likelihoods
    l0 = binom.logpmf(x, n, p0).sum() 
    l1 = binom.logpmf(x,n, p1).sum()
    # calculate the pvalue (sf=1-csf(), 1df)
    return chi2.sf(2*(l1-l0),1),(*p1,p0)

#from https://stackoverflow.com/questions/54505173/finding-alpha-and-beta-of-beta-binomial-distribution-with-scipy-optimize-and-log
def loglike_betabinom(params, k,n):
    a, b = params[0], params[1]
    logpdf = gammaln(n+1) + gammaln(k+a) + gammaln(n-k+b) + gammaln(a+b) - \
     (gammaln(k+1) + gammaln(n-k+1) + gammaln(a) + gammaln(b) + gammaln(n+a+b))
    return -np.sum(logpdf) 

def loglike_betabinom_alt(params, k,n):
    #alternative parametrization with mean, var
    a=params[0]/ params[1]-params[0]
    b=(( params[1]-1)*(params[0]-1))/ params[1]
    logpdf = gammaln(n+1) + gammaln(k+a) + gammaln(n-k+b) + gammaln(a+b) - \
     (gammaln(k+1) + gammaln(n-k+1) + gammaln(a) + gammaln(b) + gammaln(n+a+b))
    return -np.sum(logpdf) 


def betabinom_lr_test(x,n):
    # likelihood ratio test
    # x,n should be 2 sets (the two groups)
    # x is betabinomial(n,a,b), eg binomial distribution, where p follows beta ditribution with parameters a,b>0
    # mean m=a/(a+b) overdispersion g=1/(a+b+1) --> a=m/g-m and b=((g-1)*(m-1))/g
    # principle: log likelihood ratio of M0/M1 is chi2 distributed
    
    params=list()
    for xi,ni in itertools.chain(zip(x,n),((np.concatenate(x),np.concatenate(n)),)):
        xi, ni=xi[ni>0], ni[ni>0] #avoid div by 0
        #find good initialization parameters for a and b
        #prob=np.array([xii/nii for xii,nii in zip(xi,ni) if nii>0])
        prob=xi/ni
        minit=prob.mean()
        ginit=max(prob.std(),1e-6)
        init_params = [minit/ginit-minit, ((ginit-1)*(minit-1))/ginit]
        #find ml estimates for a and b
        mle = minimize(loglike_betabinom, x0=init_params,bounds=((1e-6,None),(1e-6,None)),  args=(xi,ni),options={'maxiter': 250})
        params.append(mle.x)                
    # calculate the log likelihoods
    try:
        l0 = betabinom.logpmf(np.concatenate(x), np.concatenate(n), *params[2]).sum() 
        l1 = betabinom.logpmf(x[0],n[0], *params[0]).sum()+betabinom.logpmf(x[1],n[1], *params[1]).sum()
    except ValueError:
        log.error(f'betabinom error: x={x}\nn={n}\nparams={params}')
        raise
    params_alt=[(a/(a+b), 1/(a+b+1)) for a,b in params]
    return chi2.sf(2*(l1-l0),2), params_alt #note that we need two degrees of freedom here as h0 hsa two parameters, h1 has 4

def altsplice_test(transcriptome,groups, min_cov=10, min_n=10, min_sa=.5, test=betabinom_lr_test,padj_method='fdr_bh'):
    #multitest_default={}
    sa_idx={sa:idx[0] for sa,idx in transcriptome.get_sample_idx().items()}
    assert len(groups) == 2 , "length of groups should be 2, but found %i" % len(groups)
    if isinstance(groups, dict):
        groupnames=list(groups)
        groups=list(groups.values())
    else:
        groupnames=['grp1','grp2']

    notfound=[sa for grp in groups for sa in grp if sa not in sa_idx]
    if notfound:
        raise ValueError(f"Cannot find the following samples: {notfound}")    
    grp_idx=[[sa_idx[sa] for sa in grp] for grp in groups]
    if min_sa<1:
        min_sa*=max(len(gr) for gr in groups)
    res=[]
    for g in tqdm(transcriptome):
        for junction_cov,total_cov,start,end in g.splice_graph.get_splice_coverage():
            x=[junction_cov[grp] for grp in grp_id]
            n=[total_cov[grp] for grp in grp_id]
            if any((ni>=min_n).sum()<min_sa for ni in n):
                continue
            x_sum=sum(xi.sum() for xi in x)     
            n_sum=sum(ni.sum() for ni in n)       
            if x_sum < min_cov or n_sum-x_sum < min_cov:
                continue
            pval, params=test(x,n)
            res.append(tuple(itertools.chain((g.name,g.id,g.chrom, start, end,pval),params ,
                (val for lists in zip(x,n) for pair in zip(*lists) for val in pair ))))
    df=pd.DataFrame(res, columns= (['gene','gene_id','chrom', 'start', 'end','pvalue']+ 
            ['prop_'+gn for gn in groupnames+['total']]+  
            [f'{w}_{sa}_{gn}' for gn,grp in zip(groupnames, groups) for sa in grp for w in ['cov', 'span_cov'] ]))
    mask = np.isfinite(df['pvalue'])
    padj = np.empty(mask.shape)
    padj.fill(np.nan) 
    padj[mask] = multi.multipletests(df.loc[mask,'pvalue'],method='fdr_bh')[1]
    df.insert(5,'padj',padj)
    return df

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

# summary tables (can be used as input to plot_bar / plot_dist)
def altsplice_stats(transcriptome, coverage=True,groups=None,  min_coverage=1, isoseq_filter={}):    #todo: filter like in make table
    try:
        runs=transcriptome.runs
    except KeyError:
        runs=['transcriptome']

    weights=dict()
    #if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    for _,_,tr in transcriptome.iter_transcripts(**isoseq_filter):
        w=tr['coverage'] 
        if groups is not None:
            w=[sum(w[gi] for gi in g) for g in groups.values()]
        if not coverage:
            w=[1 if wi>=min_coverage else 0 for wi in w]
        if tr['annotation'] is None:
            weights['novel/unknown']=weights.get('novel/unknown',np.zeros(len(w)))+w
        else:
            for stype in tr['annotation']['as']:
                weights[stype]=weights.get(stype,np.zeros(len(w)))+w
        weights['total']=weights.get('total',np.zeros(len(w)))+w

    df=pd.DataFrame(weights, index=runs if groups is None else groups.keys()).T
    #if groups is not None:
    #    df=pd.DataFrame({grn:df[grp].sum(1) for grn, grp in groups.items()})
        # sum per group
    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    if coverage:
        title='Expressed Transcripts'
        ylab='fraction of reads'
    else:
        title='Different Transcripts'
        ylab='fraction of  different transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'

    return df, {'ylabel':ylab,'title':title}
    #
    
def filter_stats(transcriptome, coverage=True,groups=None,  min_coverage=1,consider=['TRUNCATION', 'A_CONTENT',  'CLIPPED_ALIGNMENT', 'RTTS']):    #todo: filter like in make table
    try:
        runs=transcriptome.runs
    except KeyError:
        runs=['transcriptome']
    weights=dict()
    #if groups is not None:
    #    gi={r:i for i,r in enumerate(runs)}
    #    groups={gn:[gi[r] for r in gr] for gn,gr in groups.items()}
    for g,trid,tr in transcriptome.iter_transcripts():
        w=tr['coverage'] 
        if groups is not None:
            w=[sum(w[gi] for gi in g) for g in groups.values()]
        if not coverage:
            w=[1 if wi>=min_coverage else 0 for wi in w]
        
        relevant_filter=[f for f in tr['filter'] if f in consider]
        if relevant_filter:
            for f in relevant_filter:
                weights[f]=weights.get(f,np.zeros(len(w)))+w
        else:
            weights['PASS']=weights.get('PASS',np.zeros(len(w)))+w
        weights['total']=weights.get('total',np.zeros(len(w)))+w
            
    df=pd.DataFrame(weights, index=runs if groups is None else groups.keys()).T

    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    ylab='fraction of reads' if coverage else 'fraction of different transcripts'
    if coverage:
        title='Expressed Transcripts'
    else:
        title='Different Transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'
    return df, {'ylabel':ylab,'title':title}


def transcript_length_hist(transcriptome=None,reference=None,groups=None,bins=50,x_range=(100,10000),coverage=True,min_coverage=1,use_alignment=False, isoseq_filter={}, reference_filter={}):
    trlen=[]
    cov=[]
    #iso_filter=[]
    for _,_,tr in transcriptome.iter_transcripts(**isoseq_filter):
        cov.append(tr['coverage'])
        trlen.append(sum(e[1]-e[0] for e in tr['exons']) if use_alignment else tr['source_len'])
    #cov=pd.DataFrame(cov, columns=transcriptome.runs)
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(trlen, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if reference is not None:
        ref_len=[sum(e[1]-e[0] for e in tr['exons']) for _,_,tr in reference.iter_transcripts(**reference_filter)]
        counts['reference']=np.histogram(ref_len, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='linear', title='transcript length',xlabel='transcript length [bp]', ylabel='density', density=True)
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def transcript_coverage_hist(transcriptome, groups=None,bins=50,x_range=(0.5,1000), isoseq_filter={}, reference_filter={}):
    # get the transcript coverage in bins for groups
    # return count dataframe and suggested default parameters for plot_distr
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(**isoseq_filter):
        cov.append(tr['coverage'])
    #cov=pd.DataFrame(cov, columns=transcriptome.runs)
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(g_cov, bins=bins)[0] for gn,g_cov in cov.items()})
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='transcript coverage',xlabel='reads per transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    #plot histogram
    # cov.mask(cov.lt(x_range[0]) | cov.gt(x_range[1])).plot.hist(ax=ax, alpha=0.5, bins=n_bins)
    # ax=counts.plot.bar() 
    # ax.plot(x, counts)
   
def transcripts_per_gene_hist(transcriptome, reference=None, groups=None,bins=50,x_range=(.5,49.5), min_coverage=1, isoseq_filter={}, reference_filter={}):
    ntr=[]
    for g in transcriptome:        
        cov=np.array([g.transcripts[trid]['coverage'] for trid in g.filter_transcripts(**isoseq_filter)]).T
        ntr.append([(c>=min_coverage).sum() for c in cov])
    #ntr=pd.DataFrame(ntr, columns=transcriptome.runs)
    ntr=pd.DataFrame(ntr)
    if groups is not None:
        ntr=pd.DataFrame({grn:ntr[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(trl, bins=bins)[0] for gn,trl in ntr.items()})   
    if reference is not None:
        ref_ntr=[g.n_transcripts for g in reference]
        counts['reference']=np.histogram(ref_ntr, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    sub=f'counting transcripts covered by >= {min_coverage} reads'
    if isoseq_filter is not None:
        if 'include' in isoseq_filter:
            sub+=f'; only the following categories:{isoseq_filter["include"]}'
        if 'remove' in isoseq_filter:
            sub+=f'; without the following categories:{isoseq_filter["remove"]}'
    params=dict(yscale='log',title='transcript per gene\n'+sub,xlabel='transcript per gene', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    
def exons_per_transcript_hist(transcriptome, reference=None, groups=None,bins=35,x_range=(0.5,69.5),coverage=True,  min_coverage=2, isoseq_filter={}, reference_filter={}):
    n_exons=[]
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(**isoseq_filter):
        cov.append(tr['coverage'])
        n_exons.append(len(tr['exons']))
    #cov=pd.DataFrame(cov, columns=transcriptome.runs)
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(n_exons, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if reference is not None:
        ref_n_exons=[len(tr['exons']) for g in reference for tr in g.transcripts.values()]
        counts['reference']=np.histogram(ref_n_exons, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='exons per transcript',xlabel='number of exons per transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def downstream_a_hist(transcriptome, reference=None, groups=None,bins=30,x_range=(0,1),coverage=True,  min_coverage=2, isoseq_filter={}, reference_filter={}):
    acontent=[]
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(**isoseq_filter):
        cov.append(tr['coverage'])
        acontent.append(tr['downstream_A_content'])
    #cov=pd.DataFrame(cov, columns=transcriptome.runs)
    cov=pd.DataFrame(cov)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(x_range[0],x_range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(acontent, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if reference is not None:
        try:
            ref_acontent=[tr['downstream_A_content'] for g in reference for tr in g.transcripts.values()]
        except KeyError:
            log.warn('No A content for reference found. Run add_biases for reference first')
        else:
            counts['reference']=np.histogram(ref_acontent, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict( title='downstream genomic A content',xlabel='fraction of A downstream the transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params







