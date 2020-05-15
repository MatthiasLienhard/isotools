import isotools
import argparse
import itertools
import datetime
import pickle
import matplotlib.pyplot as plt
import argparse
import numpy as np
import pandas as pd
from  isotools.transcriptome import Transcriptome, log
from isotools.stats import *

def load_reference(args):
    if args.anno is None:
        return None
    ref_fn=args.anno
    if args.pickle:
        try:
            return Transcriptome(ref_fn+'.isotools.pkl', chromosomes=args.chrom)    
        except FileNotFoundError:
            print('no pickled file found', flush=True)
    ref=Transcriptome(ref_fn, chromosomes=args.chrom)
    if args.pickle:
        ref.save(ref_fn+'.isotools.pkl')
    
    return ref

def load_isoseq(args, reference):
    if args.pickle:
        try:
            return Transcriptome(args.bam+'_isotools.pkl', chromosomes=args.chrom)    
        except FileNotFoundError:
            print('no pickled file found', flush=True)
    isoseq=Transcriptome(args.bam, chromosomes=args.chrom)
    isoseq.annotate(reference) 
    #this looks up truncations, junction types, direct repeats at junctions and downstream genomic a content and populates the transcripts biases field
    isoseq.add_biases(args.genome)
    if args.pickle:
        isoseq.save(args.bam+'_isotools.pkl')
    return isoseq

def filter_plots(isoseq, groups, out_stem):
    log.info('filter statistics plots')
    f_stats=[]
    f_stats.append(filter_stats(isoseq, groups=groups))
    f_stats.append(filter_stats(isoseq, groups=groups, coverage=False))
    f_stats.append(filter_stats(isoseq, groups=groups, coverage=False,min_coverage=50))
    f_stats.append(filter_stats(isoseq, groups=groups, coverage=False,min_coverage=100))
    plt.rcParams["figure.figsize"] = (15,15)
    f, ax = plt.subplots(2, 2)
    ax=ax.flatten()
    for i,fs in enumerate(f_stats):
        plot_bar(fs[0],ax=ax[i],**fs[1])  
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(out_stem+'_filter_stats.png') 

def transcript_plots(isoseq,reference, groups, out_stem ):
    log.info('transcript statistics plots')
    tr_stats=[
        transcript_coverage_hist(isoseq,  groups=groups),
        transcript_length_hist(isoseq, reference=reference, groups=groups,reference_filter=dict(include=['HIGH_SUPPORT'])),
        transcripts_per_gene_hist(isoseq, reference=reference, groups=groups),
        exons_per_transcript_hist(isoseq, reference=reference, groups=groups),
        downstream_a_hist(isoseq, reference=reference, groups=groups, isoseq_filter=dict(remove=['REFERENCE', 'MULTIEXON'])),
        downstream_a_hist(isoseq, reference=reference, groups=groups, isoseq_filter=dict( remove=['NOVEL_GENE', 'UNSPLICED']))]
    tr_stats[4][1]['title']+='\nnovel single exon genes'
    tr_stats[5][1]['title']+='\nmultiexon reference genes'
    plt.rcParams["figure.figsize"] = (20,15)

    f, ax = plt.subplots(3, 2)
    ax=ax.flatten()
    for i,ts in enumerate(tr_stats):
        plot_distr(ts[0],smooth=3,ax=ax[i],**ts[1])  
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(out_stem+'_transcript_stats.png') 

def altsplice_plots(isoseq, groups, out_stem):
    altsplice = [
        altsplice_stats(isoseq, groups=groups),
        altsplice_stats(isoseq, groups=groups,coverage=False),
        altsplice_stats(isoseq, groups=groups, coverage=False,min_coverage=100)]
    plt.rcParams["figure.figsize"] = (20,15)

    f, ax = plt.subplots(3,1)
    for i,(as_counts, as_params) in enumerate(altsplice):
        plot_bar(as_counts,ax=ax[i],drop_categories=['splice_identical'], **as_params)
    plt.tight_layout(rect=[0, 0, 1, .95])
    plt.savefig(args.out+'_altsplice.png' )

def plot_diffsplice(isoseq, reference, diffsplice,gr,n,out):
    isoseq.make_index(use_name=True)
    reference.make_index(use_name=True)

    run_num={r:i for i,r in enumerate(isoseq.infos['runs'])}
    de_tab=diffsplice.head(n)
    colors=['green', 'green']
    jparams=dict(low_cov_junctions={'color':'gainsboro','lwd':1,'draw_label':False} , 
                high_cov_junctions={'color':'green','lwd':2,'draw_label':True}, 
                interest_junctions={'color':'purple','lwd':3,'draw_label':True})
    plt.rcParams["figure.figsize"] = (20,15)
    for gene_name in de_tab['gene'].unique():
        print(gene_name)
        g=isoseq[gene_name]
        f, ax = plt.subplots(3, 1)
        _=gene_track(reference[gene_name], ax=ax[0])
        ax[0].set_xlim((g.start -100, g.end +100))
        #ax[0].set_title(f'{g.name} {g.chrom}:{g.start}-{g.end}')
        #isoseq
        joi=[tuple(p) for p in de_tab.loc[de_tab['gene']==gene_name][['start','end']].values]
        for i,sn in enumerate(gr):
            _=sashimi_plot(g, ax=ax[i+1], title='isoseq '+sn, group=[run_num[r] for r in gr[sn]], text_width=150, arc_type='both', exon_color=colors[i],junctions_of_interest=joi, **jparams)
        #illumina
        #for i,sn in enumerate(gr):
        #    _=sashimi_plot_bam(g,bam_fn.format(short[sn]), ax=ax[i+3], title='illumina day1'+sn, text_width=150, exon_color=colors[i], chrom_map=chr_map,junctions_of_interest=joi, **jparams)
        
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
    parser.add_argument('--pickle', help='pickle/unpickle intermediate results', action='store_true')
    parser.add_argument('--qc_plots', help='make qc plots', action='store_true')
    parser.add_argument('--altsplice_plots', help='plot alternative splicing', action='store_true')
    parser.add_argument('--transcript_table', help='make transcript_table', action='store_true')
    parser.add_argument('--gtf_out', help='make filtered gtf', action='store_true')
    parser.add_argument('--diff', metavar='<group1/group2>',nargs='*' , help='perform differential splicing analysis')
    parser.add_argument('--chrom', nargs='*' , help='list of chromosomes to be considered')
    parser.add_argument('--diff_plots', metavar='<n>', type=int,help='make sashimi plots for <n> top differential genes')
    
    args = parser.parse_args()
    print(args)

    if args.samples:
        samples=pd.read_csv(args.samples)
        groups=dict(samples.groupby('group')['run'].apply(list))
    else: 
        groups=None
        
    
    print(f'sample group definition: {groups}') 

    reference=load_reference(args)
    isoseq=load_isoseq(args, reference=reference)
    isoseq.add_filter()
    reference.add_filter(transcript_filter={'HIGH_SUPPORT':'transcript_support_level=="1"', 'PROTEIN_CODING':'transcript_type=="protein_coding"'})

    if args.transcript_table:
        df=isoseq.transcript_table(extra_columns=['length','n_exons','exon_starts','exon_ends', 'coverage','source_len','alt_splice','filter'])
        log.info(f'writing transcript table to {args.out}_transcripts.csv')
        df.to_csv(args.out+'_transcripts.csv')

    if args.gtf_out:
        isoseq.write_gtf(args.out+'_transcripts.gtf', use_gene_name=True, remove={'A_CONTENT','RTTS','CLIPPED_ALIGNMENT'})  

    if args.qc_plots:
        filter_plots(isoseq,groups, args.out )
        transcript_plots(isoseq, reference, groups, args.out)


    if args.altsplice_plots:
        altsplice_plots(isoseq, groups, args.out)
    
    if args.diff is not None:
        for diff_cmp in args.diff:
            gr=diff_cmp.split('/')
            if len(gr) != 2:
                print('--diff argument format error: provide two groups seperated by "/" -- skipping')
                continue
            if not all(g in groups for g in gr):
                print('--diff argument format error: groupnames not found in sample table -- skipping')
                continue
            gr={gn:groups[gn] for gn in gr}
            print(f'{" vs ".join(gr)}: {" vs ".join(str(len(grp)) for grp in gr.values())} samples', flush=True)
            res=isotools.stats.altsplice_test(isoseq, list(gr.values())).sort_values('pvalue')
            sig=res.padj<.01
            print(f'{sum(sig)} differential splice sites in {len(res.loc[sig,"gene"].unique())} genes for {" vs ".join(gr)}', flush=True)
            res.to_csv(f'{args.out}_diff_{"_".join(gr)}.csv')
            if args.diff_plots is not None:
                plot_diffsplice(isoseq, reference, res,gr, args.diff_plots, args.out)
