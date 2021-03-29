'''This module is supposed to as command line script. Its outdated and needs to be adapted to the current api'''

import isotools
import argparse
import itertools
import datetime
import pickle
import matplotlib.pyplot as plt
import argparse
import numpy as np
import pandas as pd
from  isotools.transcriptome import Transcriptome, Gene
import isotools.plots
import logging
log=logging.getLogger('run_isotools')
log.setLevel(logging.DEBUG)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s: %(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.handlers=[]
log.addHandler(log_stream)



def load_reference(args):
    if args.anno is None:
        return None
    ref_fn=args.anno
    if args.pickle and not args.force_reload:
        try:
            return Transcriptome(ref_fn+'.isotools.pkl', chromosomes=args.chrom)    
        except FileNotFoundError:
            log.info('no pickled file found')
    log.info('importing reference')
    ref=Transcriptome(ref_fn, chromosomes=args.chrom)
    if args.pickle:
        log.info('saving reference as pickle file')
        ref.save(ref_fn+'.isotools.pkl')
    
    return ref

def load_isoseq(args, reference):
    if args.pickle and not args.force_reload:
        try:
            isoseq=Transcriptome(args.bam+'_isotools.pkl', chromosomes=args.chrom)    
            if args.samples:
                sample_tab=pd.read_csv(args.samples)
                runs=sample_tab['run'].tolist()
                try:
                    idx=[runs.index(r) for r in isoseq.infos['sample_table'].run]
                except ValueError:
                    raise ValueError( f'sample table does not contain infos for run {set(isoseq.infos["sample_table"].run)-set(runs)}')
                for col in sample_tab.drop('run',axis='columns').columns:
                    isoseq.infos['sample_table'][col]=sample_tab.loc[idx,col]
            return isoseq
        except FileNotFoundError:
            log.info('no pickled file found')
    log.info('importing transcripts from bam file')
    sample_tab=pd.read_csv(args.samples)
    isoseq=Transcriptome(args.bam, chromosomes=args.chrom, sample_table=sample_tab) #todo: genome filename for mutation analysis??
    log.info('annotating transcripts')
    isoseq.annotate(reference, fuzzy_junction=9) 
    fuzzy=[tr['fuzzy_junction'] for g,trid, tr in isoseq.iter_transcripts() if 'fuzzy_junction' in tr]
    log.info(f'Fixed fuzzy junctions for {len(fuzzy)} transcripts')
    #this looks up truncations, junction types, diect repeats at junctions and downstream genomic a content and populates the transcripts biases field
    log.info('adding bias information')
    isoseq.add_biases(args.genome)
    if args.pickle:
        log.info('saving transcripts as pickle file')    
        isoseq.save(args.bam+'_isotools.pkl')
    return isoseq

def filter_plots(isoseq, groups, out_stem):
    log.info('filter statistics plots')
    f_stats=[]
    f_stats.append(isoseq.filter_stats(isoseq, groups=groups))
    f_stats.append(isoseq.filter_stats(isoseq, groups=groups, coverage=False))
    f_stats.append(isoseq.filter_stats(isoseq, groups=groups, coverage=False,min_coverage=50))
    f_stats.append(isoseq.filter_stats(isoseq, groups=groups, coverage=False,min_coverage=100))
    plt.rcParams["figure.figsize"] = (12+len(groups),15)
    f, ax = plt.subplots(2, 2)
    ax=ax.flatten()
    for i,fs in enumerate(f_stats):
        isotools.plots.plot_bar(fs[0],ax=ax[i],**fs[1])  
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(out_stem+'_filter_stats.png') 

def transcript_plots(isoseq,reference, groups, out_stem ):
    log.info('transcript statistics plots')
    tr_stats=[
        isoseq.transcript_coverage_hist(isoseq,  groups=groups),
        isoseq.transcript_length_hist(isoseq, reference=reference, groups=groups,reference_filter=dict(include=['HIGH_SUPPORT'])),
        isoseq.transcripts_per_gene_hist(isoseq, reference=reference, groups=groups),
        isoseq.exons_per_transcript_hist(isoseq, reference=reference, groups=groups),
        isoseq.downstream_a_hist(isoseq, reference=reference, groups=groups, isoseq_filter=dict(remove=['REFERENCE', 'MULTIEXON'])),
        isoseq.downstream_a_hist(isoseq, reference=reference, groups=groups, isoseq_filter=dict( remove=['NOVEL_GENE', 'UNSPLICED']))]
    tr_stats[4][1]['title']+='\nnovel single exon genes'
    tr_stats[5][1]['title']+='\nmultiexon reference genes'
    plt.rcParams["figure.figsize"] = (20,15)

    f, ax = plt.subplots(3, 2)
    ax=ax.flatten()
    for i,ts in enumerate(tr_stats):
        isotools.plots.plot_distr(ts[0],smooth=3,ax=ax[i],**ts[1])  
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(out_stem+'_transcript_stats.png') 

def altsplice_plots(isoseq, groups, out_stem):
    log.info('alternative splicing statistics plots')
    altsplice = [
        isoseq.altsplice_stats(isoseq, groups=groups),
        isoseq.altsplice_stats(isoseq, groups=groups,coverage=False),
        isoseq.altsplice_stats(isoseq, groups=groups, coverage=False,min_coverage=100)]
    plt.rcParams["figure.figsize"] = (15+2*len(groups),15)

    f, ax = plt.subplots(3,1)
    for i,(as_counts, as_params) in enumerate(altsplice):
        isotools.plots.plot_bar(as_counts,ax=ax[i],drop_categories=['splice_identical'], **as_params)
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(out_stem+'_altsplice.png' )

def altsplice_examples(isoseq, n, ignore=['splice_identical']):#return the top n covered genes for each category
    examples={}
    for g in isoseq:
        if g.chrom[:3] != 'chr':
            continue
        total_cov=sum(sum(t['coverage']) for t in g.transcripts.values())
        for trid,tr in g.transcripts.items():
            cov=sum(tr['coverage'])
            try:
                s_type=tr['annotation']['as']
            except TypeError:
                s_type=['novel/unknown']
                ref_id=None
            else:
                ref_id=tr['annotation']['ref_gene_id']
            score= cov*cov/total_cov
            for cat in s_type:
                #if score > examples.get(cat,[0])[0]:
                    examples.setdefault(cat,[]).append((score, g.name, trid,ref_id,cov,total_cov))
                

    examples={k:sorted(v,key=lambda x:-x[0] ) for k,v in examples.items()}
    return{k:v[:n] for k,v in examples.items() if k not in ignore}

def plot_altsplice_examples(isoseq,groups,illu_groups,examples, out):
    nplots=len(groups)+1
    sample_idx={r:i for i,r in enumerate(isoseq.infos['sample_table'].name)}
    if illu_groups:
        illu_sample_idx={r:i for i,r in enumerate(isoseq.infos['illumina_fn'])}
        if any(gn in illu_groups for gn in groups):
            illu_groups={gn:illu_groups[gn] for gn in groups if gn in illu_groups}
        nplots += len(illu_groups) #illumina is a dict with bam filenames
    
    plt.rcParams["figure.figsize"] = (20,5*nplots)
    
    jparams=dict(low_cov_junctions={'color':'gainsboro','lwd':1,'draw_label':False} , 
            high_cov_junctions={'color':'green','lwd':2,'draw_label':True}, 
            interest_junctions={'color':'purple','lwd':3,'draw_label':True})
    exon_color='green'
    
    for cat, best_list in examples.items():        
        log.debug(cat+str(best_list))
        for i,(score, gene_name,trid, ref_id,cov,total_cov) in enumerate(best_list):
            if ref_id is not None:
                g=isoseq[ref_id]
            else:
                g=isoseq[gene_name] # novel genes (name==id)
            try:
                info=g.transcripts[trid]["annotation"]['as'][cat]
            except TypeError:
                info=list()
            log.info(f'{i+1}. best example for {cat}: {gene_name} {trid} {info}, coverage={cov} ({cov/total_cov:%})')
            f, ax = plt.subplots(nplots, 1)
            try:
                _=isotools.stats.gene_track(reference[ref_id], ax=ax[0])
            except KeyError:
                _=isotools.stats.gene_track(Gene(g.start, g.end, dict(ID=g.id, chr=g.chrom, strand=g.strand)), ax=ax[0])
            ax[0].set_xlim(g.start-100,g.end+100)
            #isoseq
            joi=[]#set joi
            offset=1
            
            if info:
                junctions=[]
                if cat=='exon skipping':
                    exons=g.transcripts[trid]['exons']
                    for pos in info:
                        idx=next(i for i,e in enumerate(exons) if e[0]>pos[0])
                        junctions.append((exons[idx-1][1],exons[idx][0]))
                    info=junctions
                elif cat=='novel exon':
                    exons=g.transcripts[trid]['exons']
                    for i,e in enumerate(exons[1:-1]):
                        if e in info:
                            junctions.extend([(exons[i][1],e[0]), (e[1],exons[i+2][0])])
                elif cat== 'novel junction':
                    junctions=info
                for pos in junctions:
                    try:
                        if len(pos)==2 and all(isinstance(x,int) for x in pos):
                            joi.append(tuple(pos)) #if this is a junction, it gets highlighed in the plot
                    except TypeError:
                        pass
            print(f'junctions of interest: {joi}')
            for sn in groups:
                _=isotools.stats.sashimi_plot(g, ax=ax[offset], title='isoseq '+sn, group=groups[sn], text_width=150, arc_type='both', exon_color=exon_color,junctions_of_interest=joi, **jparams)
                offset+=1
            #illumina
            for sn in illu_groups:
                #sn=fn.split('/')[-1].replace('.bam','')
                _=isotools.stats.sashimi_plot_bam(g, ax=ax[offset], title='illumina '+sn, group=[illu_sample_idx[r] for r in illu_groups[sn]],text_width=150, exon_color=exon_color, junctions_of_interest=joi, **jparams)
                offset+=1
            # Temporary hack: for novel genes, g.start and end is not correct!!!
            start, end=g.segment_graph[0][0],g.segment_graph[-1][1]           
            for a in ax:
                a.set_xlim((start -100, end +100))
            f.tight_layout()
            stem=f'{out}_altsplice_{cat.replace(" ","_").replace("/","_")}_{g.name}'
            plt.savefig(f'{stem}_sashimi.png')
            #zoom
            if info:
                for pos in info:
                    if isinstance(pos,int):
                        start, end=pos, pos
                    elif len(pos)==2 and all(isinstance(x,int) for x in pos):
                        if pos[1]<pos[0]:
                            start, end=sorted([pos[0], pos[0]+pos[1]])
                        else:
                            start, end=pos
                    else:
                        continue
                    for a in ax:
                        a.set_xlim((start -100, end +100))
                    ax[0].set_title(f'{g.name} {g.chrom}:{start}-{end} {cat} (cov={cov})')

        
                    
                    plt.savefig(f'{stem}_zoom_{start}_{end}_sashimi.png')
            plt.close()





def plot_diffsplice(isoseq, reference, de_tab,gr,illu_gr,out):

    nplots=len(gr)+1
    sample_idx={r:i for i,r in enumerate(isoseq.infos['sample_table'].name)}
    if illu_gr:
        illu_sample_idx={r:i for i,r in enumerate(isoseq.infos['illumina_fn'])}
        nplots += len(illu_gr) 

    
    
    exon_color='green'
    jparams=dict(low_cov_junctions={'color':'gainsboro','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':2,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':3,'draw_label':True})
    plt.rcParams["figure.figsize"] = (20,5*nplots)
    for gene_id in de_tab['gene_id'].unique():
        g=isoseq[gene_id]
        log.info(f'sashimi plot for differentially spliced gene {g.name}')        
        f, ax = plt.subplots(nplots, 1)
        try:
            _=isotools.stats.gene_track(reference[gene_id], ax=ax[0])    
        except KeyError:
            ax[0].set(frame_on=False)
            ax[0].set_title(f'Novel Gene {g.name} {g.chrom}:{g.start:,}-{g.end:,}')
        ax[0].set_xlim((g.start -100, g.end +100))
        #ax[0].set_title(f'{g.name} {g.chrom}:{g.start}-{g.end}')
        #isoseq
        joi=[tuple(p) for p in de_tab.loc[de_tab['gene_id']==gene_id][['start','end']].values]
        offset=1
        for sn in gr:
            _=isotools.stats.sashimi_plot(g, ax=ax[offset], title='isoseq '+sn, group=gr[sn], text_width=150, arc_type='both', exon_color=exon_color,junctions_of_interest=joi, **jparams)
            offset+=1
        #illumina
        for sn in illu_gr:
            #sn=fn.split('/')[-1].replace('.bam','')
            _=isotools.stats.sashimi_plot_bam(g, ax=ax[offset], title='illumina '+ sn, group=[illu_sample_idx[r] for r in illu_gr[sn]], text_width=150, exon_color=exon_color,junctions_of_interest=joi, **jparams)
            offset+=1
        # Temporary hack: for novel genes, g.start and end is not correct!!!
        start, end=g.segment_graph[0][0],g.segment_graph[-1][1]           
        for a in ax:
            a.set_xlim((start -100, end +100))
        f.tight_layout()
        plt.savefig(f'{out}_diff_{"_".join(gr)}_{g.name}_sashimi.png')
        #zoom
        for i,row in de_tab.loc[de_tab.gene==g.name].iterrows():
            if row.start>g.start and row.end<g.end:
                for a in ax:
                    a.set_xlim((row.start -1000, row.end +1000))
                ax[0].set_title(f'{g.name} {g.chrom}:{row.start}-{row.end}')
                plt.savefig(f'{out}_diff_{"_".join(gr)}_{g.name}_zoom_{row.start}_{row.end}_sashimi.png')
        plt.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='isotools',description='process isoseq bams with isotool')
    parser.add_argument('--bam', metavar='<file.bam>', help='specify isoseq aligned bam', required=True)
    parser.add_argument('--anno', metavar='<file.gtf/gff/gff3[.gz]>', help='specify reference annotation', required=True)
    parser.add_argument('--genome', metavar='<file.fasta>', help='specify reference genome file')    
    parser.add_argument('--out', metavar='</output/directory/prefix>',default='./isotools', help='specify output path and prefix')
    parser.add_argument('--samples',metavar='<samples.csv>', help='specify csv with sample / group information')
    parser.add_argument('--illu_samples',metavar='<samples.csv>', help='specify csv with illumina sample / group information')    
    parser.add_argument('--group_by',metavar='<column name>', help='specify column used for grouping the samples', default='name')
    parser.add_argument('--pickle', help='pickle/unpickle intermediate results', action='store_true')
    parser.add_argument('--qc_plots', help='make qc plots', action='store_true')
    parser.add_argument('--altsplice_stats', help='alternative splicing barplots', action='store_true')
    parser.add_argument('--transcript_table', help='make transcript_table', action='store_true')
    parser.add_argument('--gtf_out', help='make filtered gtf', action='store_true')
    parser.add_argument('--diff', metavar='<group1/group2>',nargs='*' , help='perform differential splicing analysis')
    parser.add_argument('--chrom', nargs='*' , help='list of chromosomes to be considered')
    parser.add_argument('--diff_plots', metavar='<n>', type=int,help='make sashimi plots for <n> top differential genes')
    parser.add_argument('--altsplice_plots', metavar='<n>', type=int,help='make sashimi plots for <n> top covered alternative spliced genes for each category')
    parser.add_argument('--force_reload',help='reload transcriptomes, even in presence of pickled files', action='store_true' )

    #parser.add_argument("--illumina", nargs='*',type=lambda kv: kv.split('=',2))

    
    args = parser.parse_args()

    log.debug(f'arguments: {args}')
    #illumina=dict(args.illumina) if args.illumina  else {}


    reference=load_reference(args)
    isoseq=load_isoseq(args, reference=reference)
    log.info('adding filter for isoseq')
    isoseq.add_filter()
    log.info('adding filter for reference')
    reference.add_filter(transcript_filter={'HIGH_SUPPORT':'transcript_support_level=="1"', 'PROTEIN_CODING':'transcript_type=="protein_coding"'})
    isoseq.make_index()
    reference.make_index()

    
    isoseq.fusion_table().sort_values('total_cov', ascending=False).to_csv(args.out+'_fusion.csv')
    if args.illu_samples:
        illu_samples=pd.read_csv(args.illu_samples)
        if 'illumina_fn' not in isoseq.infos or not all(sn in isoseq.infos['illumina_fn'] for sn in illu_samples.name):
            isoseq.add_illumina_coverage(dict(zip(illu_samples['name'], illu_samples['file_name'])))
            log.info('adding/overwriting illumina coverage')
            new_illu=True
        else:
            new_illu=False
        illu_groups={}
        for cat in ('name', 'group', args.group_by):
            if cat in illu_samples:
                illu_groups.update(dict(illu_samples.groupby(cat)['name'].apply(list)))
        illu_num={sn:i for i,sn in enumerate(isoseq.infos['illumina_fn'])}
    else: 
        new_illu=False
        illu_groups=illu_num={}
        
    groups=isoseq.get_sample_idx(args.group_by)
    extended_groups={}
    for d in (isoseq.get_sample_idx(col) for col in ('name', 'group', args.group_by) if col in isoseq.infos['sample_table'].columns):
        extended_groups.update(d)

    log.debug(f'sample group definition: {groups}') 
    if args.illu_samples:
        log.debug(f'illumina sample group definition: {illu_groups}\n{isoseq.infos["illumina_fn"]}\n{illu_num}') 
    
    if args.transcript_table:
        log.info(f'writing transcript table to {args.out}_transcripts.csv')
        df=isoseq.transcript_table(extra_columns=['length','n_exons','exon_starts','exon_ends', 'grop_coverage','source_len','annotation','filter'])
        df.to_csv(args.out+'_transcripts.csv')

    if args.gtf_out:
        log.info(f'writing gtf to {args.out}_transcripts.gtf')
        isoseq.write_gtf(args.out+'_transcripts.gtf', use_gene_name=True, remove={'A_CONTENT','RTTS','CLIPPED_ALIGNMENT'})  

    if args.qc_plots:
        filter_plots(isoseq,groups, args.out )
        transcript_plots(isoseq, reference,groups,  args.out)


    if args.altsplice_stats:
        altsplice_plots(isoseq, groups, args.out)
    
    if args.altsplice_plots:
        examples=altsplice_examples(isoseq, args.altsplice_plots)    
        #isoseq,reference,groups,illu_groups,examples, out    
        plot_altsplice_examples(isoseq, reference,groups,  illu_groups, examples, args.out)


    if args.diff is not None:
        for diff_cmp in args.diff:
            gr=diff_cmp.split('/')
            log.debug(f'processing {gr}')
            if len(gr) != 2:
                log.warn('--diff argument format error: provide two groups seperated by "/" -- skipping')
                continue
            if not all(gn in extended_groups for gn in gr):
                log.warning(f'--diff argument format error: group names {[gn for gn in gr if gn not in extended_groups]} not found in sample table -- skipping')
                continue
            gr={gn:extended_groups[gn] for gn in gr}
            log.info(f'testing differential splicing in {" vs ".join(gr)}: {" vs ".join(str(len(grp)) for grp in gr.values())} samples')
            res=isotools.stats.altsplice_test(isoseq, list(gr.values())).sort_values('pvalue')
            sig=res.padj<.1
            log.info(f'{sum(sig)} differential splice sites in {len(res.loc[sig,"gene"].unique())} genes for {" vs ".join(gr)}')
            res.to_csv(f'{args.out}_diff_{"_".join(gr)}.csv', index=False)
            if args.diff_plots is not None:
                if all(gn in illu_groups or gn in illu_num for gn in gr):
                    illu_gr={gn:illu_groups[gn] if gn in illu_groups else [gn] for gn in gr}
                else:
                    illu_gr=illu_groups
                sig_tab=res.head(args.diff_plots)
                if illu_gr:
                    illu_cov=list()
                    for g,jstart, jend in zip(sig_tab.gene, sig_tab.start, sig_tab.end):
                        ji=(jstart, jend)
                        #g_cov=[0,0,0,0]
                        j_cov=[{},{}]
                        cov=isoseq[g].illumina_coverage
                        for gi,grp_n in enumerate(gr):  
                            if grp_n not in illu_gr:
                                j_cov[gi]='NA'          
                            for sn in illu_gr[grp_n]:
                                i=illu_num[sn]
                                for k,v in cov[i].junctions.items():
                                    j_cov[gi][k]=j_cov[gi].get(k, 0)+v
                                j_cov[gi].setdefault(ji,0)
                        illu_cov.append((j_cov[0][ji],j_cov[1][ji], max(j_cov[0].values()), max(j_cov[1].values())))
                    illu_cov={k:v for k,v in zip(['illu_cov1','illu_cov2','illu_max1','illu_max2'],zip(*illu_cov))}
                    sig_tab=sig_tab.assign(**illu_cov)
                    
                sig_tab.to_csv(f'{args.out}_diff_top_{"_".join(gr)}.csv')    
                plot_diffsplice(isoseq, reference, res.head(args.diff_plots),gr,illu_gr, args.out)
    if args.pickle and new_illu:
        log.info('saving transcripts as pickle file (including loaded illumina profiles)')    
        isoseq.save(args.bam+'_isotools.pkl')

