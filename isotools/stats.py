from math import log10
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import pandas as pd
import numpy as np
from pysam import  AlignmentFile
from scipy.stats import binom,norm, chi2
import statsmodels.stats.multitest as multi
import isotools.splice_graph
from isotools.transcriptome import log
from tqdm import tqdm

def overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False

#sashimi plots
def sashimi_plot_bam(g, bam_fn,ax=None,text_width=.02, text_height=1, title=None,group=None, exon_color='blue', junction_color='green', min_junction_width=10, min_junction_cov=1,chrom_map=None):
    delta=np.zeros(g.end-g.start)
    chrom=g.chrom
    if chrom_map is not None:
        chrom=chrom_map[chrom]
    reg=f'{chrom}:{g.start}-{g.end}'
    if title is None:
        title=g.name
    junctions={}
    with AlignmentFile(bam_fn, "rb") as align:
        for read in align.fetch(region=reg):        
            blocks=read.get_blocks() #replace by exons?
            for i, block in enumerate(blocks):
                s=max(g.start,min(g.end-1,block[0]))-g.start
                e=max(g.start,min(g.end-1,block[1]))-g.start
                delta[s]+=1
                delta[e]-=1            
                if i>0:
                    jpos=(blocks[i-1][1],block[0])
                    if jpos[1]-jpos[0]<min_junction_width:
                        continue
                    junctions.setdefault(jpos,0)
                    junctions[jpos]+=1            
    junctions={k:v for k,v in junctions.items() if not v<min_junction_cov}
    cov=np.cumsum(delta)
    if text_width<1:
        text_width=(g.end-g.start)*text_width
    if ax is None:
        fig,ax = plt.subplots(1)
        
    ax.fill_between(range(g.start, g.end), 0, np.log10(cov+.5),facecolor=exon_color )
    
    textpositions=[]
    for (x1,x2),w in junctions.items():
        if x1<=g.start or x2>g.end:
            continue
        y1=cov[x1-g.start-1]
        y2=cov[x2-g.start]
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

        bow1=patches.Arc((center, log10(y1+.5)), width=width, height=bow_height[0]*2,theta1=90, theta2=180,linewidth=1,edgecolor=junction_color)
        bow2=patches.Arc((center, log10(y2+.5)), width=width, height=bow_height[1]*2,theta1=0, theta2=90,linewidth=1,edgecolor=junction_color)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        txt=ax.text(center,log10(max(y1,y2)+.5)+min(bow_height)+text_height/3,w,horizontalalignment='center', verticalalignment='bottom',bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
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
    ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    ax.set_title(title)
    return(ax)

def sashimi_plot(g, ax=None,text_width=.02, arc_type='coverage',text_height=1,title=None,group=None,  high_cov_th=.1,junctions_of_interest=None,  exon_color='blue', 
                low_cov_junctions={'color':'grey','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':1,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':2,'draw_label':True}):
        jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
        if 'splice_graph' not in g.data:
            g.data['splice_graph']=isotools.splice_graph.SpliceGraph(g)
        sg=g.data['splice_graph']
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
        
        idx=list(range(len(sg)))
        arcs=[]
        for i,(es,ee, pre, suc) in enumerate(sg._graph):
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
                txt=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,lab,
                        horizontalalignment='center', verticalalignment='bottom',zorder=10+priority,
                        bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
            #bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

        #ax.set_yscale('log') 
        ax.set_xlim(sg[0][0]-100, sg[-1][1]+100)
        if textpositions:
            ax.set_ylim(-text_height,max(tp[1] for tp in textpositions)+2*text_height)
        else:
            ax.set_ylim(-text_height,3) #todo: adjust y axis and ticklabels to coverage
        ax.set(frame_on=False)
        ax.set_yticks([0,1,2,3])
        ax.set_yticklabels([1,10,100,1000])
        
        def mega(x, pos):
            #'The two args are the value and tick position'
            return f'{x*1e-6} MB'

        mega_form = FuncFormatter(mega)
        ax.xaxis.set_major_formatter(mega_form)

        #ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
        #ax.set_xscale(1e-6, 'linear')
        ax.set_title(title)
        return(ax)

def gene_plot(g, ax=None,title=None, draw_exon_numbers=True, color='blue'):
    if ax is None:
        fig,ax = plt.subplots(1)    
    for i,(trid, tr) in enumerate(g.transcripts.items()):
        plt.plot(tr[exons][0][1], tr[exons][-1][0], [i+.25]*2)
        #draw arrow to mark direction
        for st, end in tr['exons']:
            rect = patches.Rectangle((st,i),(end-st),i+.5,linewidth=1,edgecolor=color,facecolor=color)
            ax.add_patch(rect)  
    return ax


# differential splicing
def proportion_test(x,n):
    # Normal approximation
    #x,n should be lenght 2(the two groups)
    #tests H0: proportions are equal vs H1: proportions are different (two sided)
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    z=abs(p1[0]-p1[1])/np.sqrt(p0*(1-p0)*(1/n[0]+1/n[1]))
    return(2*norm.sf(z)) #two sided alternative
    
def lr_test(x,n):
    # likelihood ratio test
    # x,n should be length 2 (the two groups)
    # principle: log likelihood ratio of M0/M1 is chi2 distributed
    p1=[x[i]/n[i] for i in range(2)]
    p0=(x[0]+x[1])/(n[0]+n[1])
    # calculate the log likelihoods
    l0 = binom.logpmf(x, n, p0).sum() 
    l1 = binom.logpmf(x,n, p1).sum()
    # calculate the pvalue (sf=1-csf(), 1df)
    return chi2.sf(2*(l1-l0),1)

def altsplice_test(transcriptome,groups, min_cov=10, test=proportion_test,padj_method='fdr_bh'):
    multitest_default={}
    grp_idx={r:i for i,r in enumerate(transcriptome.runs)}
    grp=[[grp_idx[r] for r in g] for g in groups]
    res=[]
    for g in tqdm(transcriptome):
        for junction_cov,total_cov,start,end in g.splice_graph.get_splice_coverage():
            x=[junction_cov[g].sum() for g in grp]
            n=[total_cov[g].sum() for g in grp]            
            if sum(x) > min_cov and sum(n)-sum(x) > min_cov and all(ni>0 for ni in n):
                p=test(x,n)
                res.append((g.name,g.chrom, start, end,p,*x,*n))
    df=pd.DataFrame(res, columns=['gene','chrom', 'start', 'end','pvalue','x1','x2','n1','n2'])
    df.insert(5,'padj',multi.multipletests(df['pvalue'],method=padj_method)[1])
    return df

#plots
def plot_bar(df,ax=None, drop_categories=None, legend=True,**axparams):    
    if ax is None:
        fig, ax=  plt.subplots()
    fractions=(df/df.sum()*100)
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

def plot_distr(counts,ax=None,density=False, legend=True,fill=True,**axparams):
    #maybe add smoothing
    x=[sum(bin)/2 for bin in counts.index]
    sz=[bin[1]-bin[0] for bin in counts.index]
    if ax is None:
        fig, ax=  plt.subplots()
    if density: 
        counts=(counts/counts.sum())
    if fill:
        for gn,gc in counts.items():
            ax.fill_between(x, 0, gc/sz, label=gn, alpha=.5)
    ax.plot(x, counts.divide(sz, axis=0))
    ax.set(**axparams)
    if legend:
        ax.legend()
    return ax

# summary tables (can be used as input to plot_bar / plot_dist)
def altsplice_stats(transcriptome, coverage=True,groups=None,  min_coverage=1, region=None, include=None, remove=None):    #todo: filter like in make table
    try:
        runs=transcriptome.runs
    except KeyError:
        runs=['transcriptome']

    weights=dict()
    for _,_,tr in transcriptome.iter_transcripts(region, include=include, remove=remove):
        w=tr['coverage'] if coverage else [1 if wi>=min_coverage else 0 for wi in tr['coverage']]
        if tr['annotation'] is None:
            weights['novel/unknown']=weights.get('novel/unknown',np.zeros(len(w)))+w
        else:
            for stype in tr['annotation']['sType']:
                weights[stype]=weights.get(stype,np.zeros(len(w)))+w
                    
    df=pd.DataFrame(weights, index=runs).T
    if groups is not None:
        df=pd.DataFrame({grn:df[grp].sum(1) for grn, grp in groups.items()})
        # sum per group
    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)

    return df,  {'ylabel':'fraction of reads' if coverage else 'fraction of different transcripts'}
    #
    
def filter_stats(transcriptome, coverage=True,groups=None,  min_coverage=1,consider=['TRUNCATION', 'A_CONTENT',  'CLIPPED_ALIGNMENT', 'RTTS']):    #todo: filter like in make table
    try:
        runs=transcriptome.runs
    except KeyError:
        runs=['transcriptome']
    weights=dict()
    for g,trid,tr in transcriptome.iter_transcripts():
        w=tr['coverage'] if coverage else [1 if wi>=min_coverage else 0 for wi in tr['coverage']]
        relevant_filter=[f for f in tr['filter'] if f in consider]
        if relevant_filter:
            for f in relevant_filter:
                weights[f]=weights.get(f,np.zeros(len(w)))+w
        else:
            weights['PASS']=weights.get('PASS',np.zeros(len(w)))+w
            
    df=pd.DataFrame(weights, index=runs).T
    if groups is not None:
        df=pd.DataFrame({grn:df[grp].sum(1) for grn, grp in groups.items()})
        # sum per group
    df=df.reindex(df.mean(1).sort_values(ascending=False).index, axis=0)
    ylab='fraction of reads' if coverage else 'fraction of different transcripts'
    if coverage:
        title='Expressed Transcripts'
    else:
        title='Different Transcripts'
        if min_coverage>1:
            title+=f' > {min_coverage} reads'
    return df, {'ylabel':ylab,'title':title}
   
    if ax is None:
        fig, ax=  plt.subplots()
    fractions=(df/df.sum()*100)
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

def transcript_length_hist(transcriptome,reference=None,groups=None,bins=50,range=(100,10000),coverage=True,min_coverage=1,use_alignment=False, include=None, remove=None):
    trlen=[]
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(include=include, remove=remove):
        cov.append(tr['coverage'])
        trlen.append(sum(e[1]-e[0] for e in tr['exons']) if use_alignment else tr['source_len'])
    cov=pd.DataFrame(cov, columns=transcriptome.runs)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(range[0],range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(trlen, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if reference is not None:
        ref_len=[sum(e[1]-e[0] for e in tr['exons']) for g in reference for tr in g.transcripts.values()]
        counts['reference']=np.histogram(ref_len, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='linear', title='transcript length',xlabel='transcript length [bp]', ylabel='density', density=True)
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

def transcript_coverage_hist(transcriptome, groups=None,bins=50,range=(0.5,1000), include=None, remove=None):
    # get the transcript coverage in bins for groups
    # return count dataframe and suggested default parameters for plot_distr
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(include=include, remove=remove):
        cov.append(tr['coverage'])
    cov=pd.DataFrame(cov, columns=transcriptome.runs)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(range[0],range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(g_cov, bins=bins)[0] for gn,g_cov in cov.items()})
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='transcript coverage',xlabel='reads per transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    #plot histogram
    # cov.mask(cov.lt(range[0]) | cov.gt(range[1])).plot.hist(ax=ax, alpha=0.5, bins=n_bins)
    # ax=counts.plot.bar() 
    # ax.plot(x, counts)
   
def transcripts_per_gene_hist(transcriptome, reference=None, groups=None,bins=50,range=(.5,49.5), min_coverage=1,include=None, remove=None):
    ntr=[]
    for g in transcriptome:        
        cov=np.array([g.transcripts[trid]['coverage'] for trid in g.filter_transcripts(include, remove)]).T
        ntr.append([(c>=min_coverage).sum() for c in cov])
    ntr=pd.DataFrame(ntr, columns=transcriptome.runs)
    if groups is not None:
        ntr=pd.DataFrame({grn:ntr[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(range[0],range[1],bins)
    counts=pd.DataFrame({gn:np.histogram(trl, bins=bins)[0] for gn,trl in ntr.items()})   
    if reference is not None:
        ref_ntr=[g.n_transcripts for g in reference]
        counts['reference']=np.histogram(ref_ntr, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    sub=f'counting transcripts covered by >= {min_coverage} reads'
    if include is not None:
        sub+=f'; only the following categories:{include}'
    if include is not None:
        sub+=f'; without the following categories:{remove}'
    params=dict(yscale='log',title='transcript per gene\n'+sub,xlabel='transcript per gene', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params
    
def exons_per_transcript_hist(transcriptome, reference=None, groups=None,bins=35,range=(0.5,69.5),coverage=True,  min_coverage=2,include=None, remove=None ):
    n_exons=[]
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(include=include, remove=remove):
        cov.append(tr['coverage'])
        n_exons.append(len(tr['exons']))
    cov=pd.DataFrame(cov, columns=transcriptome.runs)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(range[0],range[1],bins)
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


def downstream_a_hist(transcriptome, reference=None, groups=None,bins=30,range=(0,1),coverage=True,  min_coverage=2,include=None, remove=None ):
    acontent=[]
    cov=[]
    for _,_,tr in transcriptome.iter_transcripts(include=include, remove=remove):
        cov.append(tr['coverage'])
        acontent.append(tr['biases']['downstream_A_content'])
    cov=pd.DataFrame(cov, columns=transcriptome.runs)
    if groups is not None:
        cov=pd.DataFrame({grn:cov[grp].sum(1) for grn, grp in groups.items()})
    if isinstance(bins,int):
        bins=np.linspace(range[0],range[1],bins)
    if not coverage:
        cov[cov<min_coverage]=0
        cov[cov>0]=1
    counts=pd.DataFrame({gn:np.histogram(acontent, weights=g_cov, bins=bins)[0] for gn,g_cov in cov.items()})   
    if reference is not None:
        try:
            ref_acontent=[tr['biases']['downstream_A_content'] for g in reference for tr in g.transcripts.values()]
        except KeyError:
            log.warn('No A content for reference found. Run add_biases for reference first')
        else:
            counts['reference']=np.histogram(ref_acontent, bins=bins)[0]
    bin_df=pd.DataFrame({'from':bins[:-1],'to':bins[1:]})
    params=dict(yscale='log', title='downstream genomic A content',xlabel='fraction of A downstream the transcript', ylabel='# transcripts')
    return pd.concat([bin_df,counts], axis=1).set_index(['from', 'to']),params

#biases plots






