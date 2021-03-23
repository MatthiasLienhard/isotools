  
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

def sashimi_figure(self,samples=None,short_read_samples=None,draw_gene_track=True,long_read_params=None, short_read_params=None ,junctions_of_interest=None, x_range=None):
    draw_gene_track=bool(draw_gene_track)
    if samples is None and short_read_samples is None:
        samples={'all':self._transcriptome.samples} #all samples grouped # pylint: disable=W0212
    if samples is None:
        samples={}
    else:
        samples=get_index(samples, self._transcriptome.samples)# pylint: disable=W0212
    if short_read_samples is None:
        short_read_samples={}
    else:
        short_read_samples=get_index(short_read_samples, self._transcriptome.infos['short_reads']['name'])    # pylint: disable=W0212
    long_read_params=extend_params(long_read_params)
    short_read_params=extend_params(short_read_params)
    if x_range is None:
        x_range=(self.start-100, self.end+100)

    f,axes=plt.subplots(len(samples)+len(short_read_samples) +draw_gene_track)
    axes=np.atleast_1d(axes) # in case there was only one subplot

    if draw_gene_track:
        self.gene_track(ax=axes[0],x_range=x_range)        

    for i,(sname,sidx) in enumerate(samples.items()):
        self.sashimi_plot(sidx,sname,axes[i+draw_gene_track],junctions_of_interest,x_range=x_range, **long_read_params)
        

    for i,(sname,sidx) in enumerate(short_read_samples.items()):
        self.sashimi_plot_short_reads(sidx,sname,axes[i+len(samples)+draw_gene_track],junctions_of_interest,x_range=x_range, **long_read_params)
        
    
    return f,axes


def sashimi_plot_short_reads(self,sidx, title, ax,junctions_of_interest, x_range,jparams, exon_color, high_cov_th,min_cov_th, text_width, arc_type,text_height):
    short_reads=[self.short_reads(idx) for idx in sidx]
    #jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    start=short_reads[0].reg[1]
    end=short_reads[0].reg[2]
    #delta=np.zeros(end-start)
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
            _=ax.text(center,log10(max(y1,y2))+min(bow_height)+text_height/3,w,horizontalalignment='center', verticalalignment='bottom',bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)).set_clip_on(True)
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

def sashimi_plot(self,sidx, title, ax,junctions_of_interest, x_range,jparams, exon_color, high_cov_th,min_cov_th, text_width, arc_type,text_height):   
    ebg=self.eb_graph
    boxes=[(node[0], node[1], self.coverage[np.ix_(sidx,[i for i in set(node[2]).union(node[3])])].sum()) for node in ebg]
    if text_width<1:
        text_width=(ebg[-1][1]-ebg[0][0])*text_width
    total_weight=self.coverage[sidx,:].sum()
    if high_cov_th<1:
        high_cov_th=high_cov_th*total_weight   
    if min_cov_th<1:
        min_cov_th=min_cov_th*total_weight   
    #idx=list(range(len(ebg)))
    arcs=[]
    for i,(_,ee, _, suc) in enumerate(ebg):
        weights={}
        for tr,next_i in suc.items():
            weights.setdefault(next_i,0)
            weights[next_i]+=self.coverage[np.ix_(sidx,[tr])].sum()
        arcs_new=[(ee,boxes[i][2],ebg[next_i][0],boxes[next_i][2],w) for next_i, w in weights.items() if ebg[next_i][0]>ee and w>0]
        if arcs_new:
            arcs.extend(arcs_new)
    if ax is None:
        _,ax = plt.subplots(1)
        
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

def gene_track(self, ax=None,title=None, reference=True, remove_transcripts=None, label_exon_numbers=True,label_transcripts=True,label_fontsize=10, color='blue',x_range=None):
    contrast='white' if np.mean(plt_col.to_rgb(color))<.5 else 'black'
    if ax is None:
        _,ax = plt.subplots(1)    
    if x_range is None:
        x_range=(self.start-100, self.end+100)
    blocked=[]
    if remove_transcripts is None:
        remove_transcripts=[]
    if reference: #select transcripts and sort by start
        transcript_list=sorted([(tr_nr,tr) for tr_nr,tr in enumerate(self.ref_transcripts) if tr_nr not in remove_transcripts],key=lambda x:x[1]['exons'][0][0]) #sort by start position
    else:
        transcript_list=sorted([(tr_nr,tr) for tr_nr,tr in enumerate(self.transcripts)     if tr_nr not in remove_transcripts],key=lambda x:x[1]['exons'][0][0]) 
    for tr_nr,tr in transcript_list:
        tr_start, tr_end=tr['exons'][0][0],tr['exons'][-1][1]
        if (tr_end < x_range[0] or tr_start > x_range[1]): #transcript does not overlap x_range
            continue
        trid='> ' if self.strand=='+' else '< ' # indicate the strand like in ensembl browser
        trid+=tr['transcript_name'] if 'transcript_name' in tr else f'{self.id}_{tr_nr}'
        
        # find next line that is not blocked
        try:
            i=next(idx for idx,last in enumerate(blocked) if last<tr['exons'][0][0] )
        except StopIteration:
            i=len(blocked)
            blocked.append(tr_end)
        else:
            blocked[i]=tr_end
        #line from TSS to PAS at 0.25
        ax.plot((tr_start, tr_end), [i+.25]*2, color=color)
        if label_transcripts:
            pos=(max(tr_start, x_range[0])+min(tr_end, x_range[1]))/2
            ax.text(pos,i-.02,trid, ha='center', va='top',fontsize=label_fontsize, clip_on=True)
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
            if label_exon_numbers and (end>x_range[0] and st<x_range[1]):
                enr=j+1 if self.strand=='+' else len(tr['exons'])-j
                pos=(max(st, x_range[0])+min(end, x_range[1]))/2
                ax.text(pos,i+.25,enr,ha='center', va='center', color=contrast, fontsize=label_fontsize, clip_on=True)    #bbox=dict(boxstyle='round', facecolor='wheat',edgecolor=None,  alpha=0.5)
        i+=1
    if title is None:
        title=self.name
    ax.set_title(title)
    ax.set(frame_on=False)   
    ax.get_yaxis().set_visible(False)
    ax.set_ylim(-.5,len(blocked))
    ax.set_xlim(*x_range)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    return ax
