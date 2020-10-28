  
import matplotlib.colors as plt_col
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np


#sashimi plots
def sashimi_plot_bam(self,ax=None,text_width=.02, text_height=1, title=None,group=None, high_cov_th=.1,chrom_map=None,junctions_of_interest=None, exon_color='blue', 
                low_cov_junctions={'color':'grey','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':1,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':2,'draw_label':True}):
    jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    delta=np.zeros(self.end-self.start)
    chrom=self.chrom
    if title is None:
        title=self.name    
    if text_width<1:
        text_width=(self.end-self.start)*text_width
    #cov=np.array([sum(illu[pos] for i, illu in enumerate(self.illumina_coverage) if groups is None or i in groups) for pos in range(self.start, self.end)])
    cov=np.zeros(self.end-self.start)
    junctions={}
    for i,illu in enumerate(self.illumina_coverage):
        if group is None or i in group:
            cov+=illu.profile
            for k,v in illu.junctions.items():            
                junctions[k]=junctions.get(k,0)+v     
                   
    total_weight=max(cov)
    if high_cov_th<1:
        high_cov_th=high_cov_th*total_weight

    if ax is None:
        fig,ax = plt.subplots(1)
        
    ax.fill_between(range(self.start, self.end), 0, np.log10(cov, where=cov>0, out=np.nan*cov),facecolor=exon_color )
    
    textpositions=[]
    for (x1,x2),w in junctions.items():
        if x1<=self.start or x2>=self.end:
            #todo: this seems to happen at some chrMT genes?
            logging.debug(f'attempt to plot junction ({(x1,x2)}) outside of gene: {self.__str__()}')
            continue
        y1=cov[x1-self.start-1]+.5
        y2=cov[x2-self.start]+.5
        center=(x1+x2)/2 
        width=x2-x1
        bow_height=text_height
        while any(label_overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
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

    
    ax.set_xlim(self.start-100, self.end+100)
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

def sashimi_plot(self, ax=None,text_width=.02, arc_type='coverage',text_height=1,
                title=None,group=None,  high_cov_th=.1,junctions_of_interest=None,  exon_color='blue', 
                low_cov_junctions={'color':'grey','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':1,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':2,'draw_label':True}):
    jparams=[low_cov_junctions,high_cov_junctions,interest_junctions]
    
    sg=self.splice_graph
    if title is None:
        title=self.name
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
        while any(label_overlap((center,log10(max(y1,y2))+bow_height), tp,text_width,text_height) for tp in textpositions):
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
    ax.set_xlim(self.start-100, self.end+100)
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

def gene_track(self, ax=None,title=None, draw_exon_numbers=True, color='blue'):
    contrast='white' if np.mean(plt_col.to_rgb(color))<.5 else 'black'
    if ax is None:
        _,ax = plt.subplots(1)    
    i=0
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
    ax.set_yticks([i+.25 for i in range(self.n_ref_transcript)])
    ax.set_yticklabels(tr['transcript_name'] for tr in self.ref_transcripts) #todo: refence transcript names
    ax.tick_params(left=False)
    ax.set_ylim(-.5,self.n_ref_transcripts+1)
    ax.set_xlim(self.start-100, self.end+100)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos=None: f'{x:,.0f}'))
    return ax
