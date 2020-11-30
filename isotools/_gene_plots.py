  
import matplotlib.colors as plt_col
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from math import log10
import logging
logger=logging.getLogger('isotools')

def _label_overlap(pos1,pos2,width, height):
    if abs(pos1[0]-pos2[0])<width and abs(pos1[1]-pos2[1])<height:
        return True
    return False


DEFAULT_JPARAMS=[   {'color':'lightgrey'  ,'lwd':1,'draw_label':False}, #low coverage junctions
                    {'color':'green' ,'lwd':1,'draw_label':True},  #high coverage junctions
                    {'color':'purple','lwd':2,'draw_label':True}]  #junctions of interest
DEFAULT_PARAMS=dict(min_cov_th=.001,high_cov_th=.05, text_width=.02, arc_type='both',text_height=1, exon_color='green')

def extend_params(params):
    if params is None:
        params=dict()
    params.setdefault('jparams',[{},{},{}])
    #jparams=[params.pop(k,jparams[i]) for i,k in enumerate(['low_cov_junctions','high_cov_junctions','interest_junctions'])]
    for i,k1 in enumerate(['low_cov_junctions','high_cov_junctions','interest_junctions']):
        params['jparams'][i]=params.pop(k1,params['jparams'][i])
        for k2,v in DEFAULT_JPARAMS[i].items():
            params['jparams'][i].setdefault(k2,v)
    for k,v in DEFAULT_PARAMS.items():
        params.setdefault(k,v)
    return params

def get_index(samples,names):
    if not samples:
        return {}
    if isinstance(names,list):
        idx={sa:i for i,sa in enumerate(names)}
    else:
        idx={sa:i for i,sa in names.items()}
    if isinstance(samples,list):
        samples=[[s] if isinstance(s,str) else s for s in samples]
        samples={','.join(sl):sl for sl in samples}
    try:
        samples={n:[idx[sa] for sa in sl] for n,sl in samples.items()}
    except KeyError:
        notfound=list({sa for sl in samples for sa in sl if sa not in idx})
        logger.error('did not find the following samples: %s',','.join(notfound))
        raise
    return samples

#sashimi plots

def sashimi_plot(self,samples=None,short_read_samples=None,gene_track=True,long_read_params=None, short_read_params=None ,junctions_of_interest=None, x_range=None):
    gene_track=bool(gene_track)
    if samples is None and short_read_samples is None:
        samples={'all':self._transcriptome.samples} #all samples grouped
    if samples is None:
        samples={}
    else:
        samples=get_index(samples, self._transcriptome.samples)
    if short_read_samples is None:
        short_read_samples={}
    else:
        short_read_samples=get_index(short_read_samples, self._transcriptome.infos['short_reads']['name'])    
    long_read_params=extend_params(long_read_params)
    short_read_params=extend_params(short_read_params)
    if x_range is None:
        x_range=(self.start-100, self.end+100)

    f,axes=plt.subplots(len(samples)+len(short_read_samples) +gene_track)

    if gene_track:
        self.gene_track(ax=axes[0],x_range=x_range)
        

    for i,(sname,sidx) in enumerate(samples.items()):
        _sashimi_plot_long_reads(self.splice_graph,sidx,sname,axes[i+gene_track],junctions_of_interest,x_range=x_range, **long_read_params)
        

    for i,(sname,sidx) in enumerate(short_read_samples.items()):
        _sashimi_plot_short_reads([self.short_reads(idx) for idx in sidx],sname,axes[i+len(samples)+gene_track],junctions_of_interest,x_range=x_range, **long_read_params)
        
    
    return f,axes


def _sashimi_plot_short_reads(short_reads,      title, ax,junctions_of_interest, x_range,jparams, exon_color, high_cov_th,min_cov_th, text_width, arc_type,text_height):
    #jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    start=short_reads[0].reg[1]
    end=short_reads[0].reg[2]
    delta=np.zeros(end-start)
    cov=np.zeros(end-start)
    junctions={}
    for sr_cov in short_reads:
        cov+=sr_cov.profile
        for k,v in sr_cov.junctions.items():            
            junctions[k]=junctions.get(k,0)+v                      
    if high_cov_th<1:
        high_cov_th*=max(cov)
    if min_cov_th<1:
        min_cov_th*=max(cov)
    #exons
    ax.fill_between(range(start, end), 0, np.log10(cov, where=cov>0, out=np.nan*cov),facecolor=exon_color )    
    #junctions
    textpositions=[]
    for (x1,x2),w in junctions.items():
        if junctions_of_interest is not None and (x1,x2) in junctions_of_interest:
            priority=2
        elif w<min_cov_th:
            continue
        elif w< high_cov_th:
            priority=0
        else:
            priority=1
        y1=cov[x1-start-1]+.5
        y2=cov[x2-start]+.5
        center=(x1+x2)/2 
        width=x2-x1
        bow_height=text_height
        while any(_label_overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
            bow_height+=text_height
        textpositions.append((center, log10(max(y1,y2))+bow_height))
        if y1<y2: 
             bow_height=(log10(y2/y1)+bow_height,bow_height)
        elif y1>y2:
             bow_height=(bow_height,bow_height+log10(y1/y2))
        else:
            bow_height=(bow_height,bow_height)
        bow1=patches.Arc((center, log10(y1)), width=width, height=bow_height[0]*2,theta1=90, theta2=180,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        bow2=patches.Arc((center, log10(y2)), width=width, height=bow_height[1]*2,theta1=0, theta2=90,linewidth=jparams[priority]['lwd'],edgecolor=jparams[priority]['color'],zorder=priority)
        ax.add_patch(bow1)
        ax.add_patch(bow2)
        if jparams[priority]['draw_label']:
            txt=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,w,horizontalalignment='center', verticalalignment='bottom',bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
        #bbox_list.append(txt.get_tightbbox(renderer = fig.canvas.renderer))

    
    ax.set_xlim(*x_range)
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

def _sashimi_plot_long_reads(splice_graph,sidx, title, ax,junctions_of_interest, x_range,jparams, exon_color, high_cov_th,min_cov_th, text_width, arc_type,text_height):   
    boxes=[(node[0], node[1], splice_graph.weights[np.ix_(sidx,[i for i in set(node[2]).union(node[3])])].sum()) for node in splice_graph]
    if text_width<1:
        text_width=(splice_graph[-1][1]-splice_graph[0][0])*text_width
    total_weight=splice_graph.weights[sidx,:].sum()
    if high_cov_th<1:
        high_cov_th=high_cov_th*total_weight   
    if min_cov_th<1:
        min_cov_th=min_cov_th*total_weight   
    #idx=list(range(len(splice_graph)))
    arcs=[]
    for i,(_,ee, _, suc) in enumerate(splice_graph._graph):
        weights={}
        for tr,next_i in suc.items():
            weights.setdefault(next_i,0)
            weights[next_i]+=splice_graph.weights[np.ix_(sidx,[tr])].sum()
        arcs_new=[(ee,boxes[i][2],splice_graph[next_i][0],boxes[next_i][2],w) for next_i, w in weights.items() if splice_graph[next_i][0]>ee and w>0]
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
        if junctions_of_interest is not None and (x1,x2) in junctions_of_interest:
            priority=2
        elif w< min_cov_th:
            continue
        elif w< high_cov_th:
            priority=0
        else:
            priority=1
        center=(x1+x2)/2 
        width=x2-x1        
        bow_height=text_height        
        while any(_label_overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
            bow_height+=text_height
        textpositions.append((center, log10(max(y1,y2))+bow_height))
        if y1<y2: 
            bow_height=(log10(y2/y1)+bow_height,bow_height)
        elif y1>y2:
            bow_height=(bow_height,bow_height+log10(y1/y2))
        else:
            bow_height=(bow_height,bow_height)
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
    if textpositions:
        ax.set_ylim(-text_height,max(tp[1] for tp in textpositions)+2*text_height)
    else:
        ax.set_ylim(-text_height,3) #todo: adjust y axis and ticklabels to coverage
    ax.set_xlim(*x_range)
    ax.set(frame_on=False)
    ax.set_yticks([0,1,2,3])
    ax.set_yticklabels([1,10,100,1000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    #ax.ticklabel_format(axis='x', style='sci',scilimits=(6,6))
    #ax.set_xscale(1e-6, 'linear')
    ax.set_title(title)

def gene_track(self, ax=None,title=None, draw_exon_numbers=True, color='blue',x_range=None):
    contrast='white' if np.mean(plt_col.to_rgb(color))<.5 else 'black'
    if ax is None:
        _,ax = plt.subplots(1)    
    i=0
    if x_range is None:
        x_range=(self.start-100, self.end+100)
    for tr in self.ref_transcripts:
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
                enr=j+1 if self.strand=='+' else len(tr['exons'])-j
                ax.text((st+end)/2,i+.25,enr,ha='center', va='center', color=contrast).set_clip_on(True)    #bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)
        i+=1
    if title is None:
        title=self.name
    ax.set_title(title)
    ax.set(frame_on=False)   
    ax.set_yticks([i+.25 for i in range(self.n_ref_transcripts)])
    ax.set_yticklabels(tr['transcript_name'] if 'transcript_name' in tr else f'transcript {i}' for i,tr in enumerate(self.ref_transcripts) ) 
    ax.tick_params(left=False)
    ax.set_ylim(-.5,self.n_ref_transcripts+1)
    ax.set_xlim(*x_range)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    return ax
