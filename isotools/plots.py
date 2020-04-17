from math import log10
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import numpy as np
from pysam import  AlignmentFile

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