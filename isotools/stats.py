from math import log10
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd
import numpy as np
from pysam import  AlignmentFile
from scipy.stats import binom,norm, chi2
import statsmodels.stats.multitest as multi
import splice_graph
from tqdm import tqdm

def overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False

def sashimi_plot_bam(g, bam_fn,ax=None,text_width=.02, text_height=1, title=None,group=None, exon_color='blue', junction_color='green', min_junction_width=10, min_junction_cov=1):
    delta=np.zeros(g.end-g.start)
    reg=f'{g.chrom}:{g.start}-{g.end}'
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

def sashimi_plot(g, ax=None,text_width=.02, arc_type='coverage',text_height=1,title=None,group=None,  exon_color='blue', junction_color='green',  min_junction_cov=1):
        if 'splice_graph' not in g.data:
            g.data['splice_graph']=splice_graph.SpliceGraph(g)
        sg=g.data['splice_graph']
        if title is None:
            title=g.name
        if group is None:
            group=list(range(sg.weights.shape[0])) #all
        
        boxes=[(node[0], node[1], sg.weights[np.ix_(group,[i for i in set(node[2]).union(node[3])])].sum()) for node in sg]
        if text_width<1:
            text_width=(sg[-1][1]-sg[0][0])*text_width
        total_weight=sg.weights[group,:].sum()
        if min_junction_cov<1:
            min_junction_cov=min_junction_cov*total_weight
        
        idx=list(range(len(sg)))
        arcs=[]
        for i,(es,ee, pre, suc) in enumerate(sg._graph):
            weights={}
            for tr,next_i in suc.items():
                weights.setdefault(next_i,0)
                weights[next_i]+=sg.weights[np.ix_(group,[tr])].sum()
            arcs_new=[(ee,boxes[i][2],sg[next_i][0],boxes[next_i][2],w) for next_i, w in weights.items() if sg[next_i][0]>ee and w>=min_junction_cov]
            if arcs_new:
                arcs.extend(arcs_new)
        if ax is None:
            fig,ax = plt.subplots(1)
            
        for st, end, h in boxes:
            if h>0:
                rect = patches.Rectangle((st,0),(end-st),log10(h),linewidth=1,edgecolor=exon_color,facecolor=exon_color)
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

            bow1=patches.Arc((center, log10(y1)), width=width, height=bow_height[0]*2,theta1=90, theta2=180,linewidth=1,edgecolor=junction_color)
            bow2=patches.Arc((center, log10(y2)), width=width, height=bow_height[1]*2,theta1=0, theta2=90,linewidth=1,edgecolor=junction_color)
            ax.add_patch(bow1)
            ax.add_patch(bow2)
            if arc_type=='coverage':
                lab=str(w)
            else:                #fraction
                lab=f'{w/total_weight:.1%}'
                if arc_type=='both':
                    lab=str(w)+' / '+lab
            txt=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,lab,horizontalalignment='center', verticalalignment='bottom',bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
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
            return '%1.3fM' % (x*1e-6)

        mega_form = FuncFormatter(mega)
        ax.xaxis.set_major_formatter(mega_form)

        #ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
        #ax.set_xscale(1e-6, 'linear')
        ax.set_title(title)
        return(ax)
        
def plot_altsplice(transcriptome,plot=True, coverage=True,groups=None,  min_coverage=1,drop_categories=None, region=None, include=None, remove=None,  ax=None):    #todo: filter like in make table
    try:
        runs=transcriptome.infos['runs']
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

    if not plot:
        return None, df
    #
    ax=plot_fractions(df=df, ax=ax, drop_categories=drop_categories)
    ax.set_ylabel('fraction of reads' if coverage else 'fraction of different transcripts')
    return ax,df

def plot_filter(transcriptome,plot=True, coverage=True,groups=None,  min_coverage=1,drop_categories=None, consider=['TRUNCATION', 'A_CONTENT',  'CLIPPED_ALIGNMENT', 'RTTS'],  ax=None):    #todo: filter like in make table
    
    try:
        runs=transcriptome.infos['runs']
    except KeyError:
        runs=['transcriptome']
    weights=dict()
    for g,trid,tr in transcriptome.iter_transcripts():
        for tr in g.transcripts.values():
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
    if not plot:
        return None, df
    #
    ax=plot_fractions(df=df, ax=ax, drop_categories=drop_categories)
    ax.set_ylabel('fraction of reads' if coverage else 'fraction of different transcripts')
    return ax,df

def plot_fractions(df,ax=None, drop_categories=None):    
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
    return ax

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
    grp_idx={r:i for i,r in enumerate(transcriptome.infos['runs'])}
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


#QC plots
def plot_transcript_len(transcriptome, plot=True, coverage=True,groups=None,min_coverage=1, include=None, remove=None):
    pass

def plot_transcript_coverage(transcriptome, plot=True, groups=None, include=None, remove=None):
    pass

def plot_transcripts_per_gene(transcriptome, plot=True, groups=None, min_coverage=1,include=None, remove=None ):
    pass

def plot_exons_per_transcript(transcriptome, plot=True, groups=None, min_coverage=1,include=None, remove=None ):
    pass

#biases plots






