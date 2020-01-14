
import os
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt

def strip(text):
    try:
        return text.strip()
    except AttributeError:
        return text

def make_int(text):
    i=text.find('(')
    if i<0: 
        i=None
    try:
        val=int(text[:i].strip('" '))
    except ValueError:
        val=None
    return val


#dir='/project/42/pacbio/hecatos_isoseq/'
def get_stats(dir, suffix,subfolder='', sumGroups=False):
    tabs=[]
    groupn=dict()
    for root, dirs, files in os.walk(dir+subfolder):
        for file in files:
            if file.endswith(suffix):
                print(file)
                group=file[:file.find('_')]
                groupn[group]=groupn.get(group,0)+1
                #tab=pd.read_csv(os.path.join(root, file))
                sampleid='{}_{}'.format(group, groupn[group])
                tabs.append(pd.read_csv(os.path.join(root, file), sep=':',
                      names=["Description", sampleid],
                      converters = {'Description' : strip,
                                    sampleid : make_int,
                                    }).set_index('Description' ))
    df=pd.concat(tabs, axis=1, sort=False)
    if sumGroups:
        colnames={gn:[t for t in df.columns if t.startswith(gn)] for gn in groupn}
        tabs=[]
        for gn, coln in colnames.items():
            tabs.append(df[coln].sum(axis=1).rename(gn))
        groupdf = pd.concat(tabs, axis=1, keys=[s.name for s in tabs])
        groupdf.loc['SMRT Cells']=[groupn[n] for n in groupdf.columns]
        return df, groupdf
    else:
        return df

def plot_altsplice(df, out_fn='altsplice.png'):    

    stype=[{'novel/unknown'} if t=='NA' else t for t in df['splice_type']]
    
    type_list=list(itertools.chain.from_iterable(stype))
    type_counts=sorted(list(zip(*np.unique(type_list, return_counts=True))), key=lambda x:x[1])
    total=len(stype)
    type_fraction=[(n,round(v/total*100,2)) for n,v in type_counts]
    D=dict(type_fraction)
    
    plt.barh(range(len(D)), list(D.values()), align='center')
    plt.yticks(range(len(D)), list(D.keys()))
    plt.savefig(out_fn,bbox_inches='tight')
    plt.close()
    for n,v in type_counts:
        print('{}:\t{}\t{}'.format(n,v, D[n]))

def plot_more(df,out_fn='plot.png', refseq_df=None):
    print('total: {}'.format(len(df)))
    stats=dict()
    stats['isoforms']=df.shape[0]
    genes=df['gene_name']

    stats['genes']=len(np.unique(genes))
    stats['refseq_isoforms']=sum(df['Splice_IOU']==1)
    novel_isoforms=df.loc[df['transcript'].isna(),'run']
    stats['novel_isoforms']=len(novel_isoforms)
    novel_genes=dict(zip(*np.unique([x[1] for x in novel_isoforms.str.split('.') if len(x)>1], return_counts=True)))
    stats['novel_genes']=len(novel_genes)
    high_isoform_genes=sorted(novel_genes.items(), key=lambda kv: kv[1], reverse=True)
    stats['high_isoform_genes']=high_isoform_genes[:10]
    for s in stats.items():
        print('{}: {}'.format(*s))
    if refseq_df is not None:
        #refseq_df=pd.read_csv(refseq_fn, sep='\t',   usecols=["name", "exonCount", "name2"], index_col='name')
        refseq_df['type']=[r[:2] for r in refseq_df.index]
        refseq_df=refseq_df.loc[refseq_df['type']=='NM']
        n_refseq_iso=np.unique(refseq_df['name2'], return_counts=True)
        n_iso=np.unique(genes, return_counts=True)
        kwargs = dict(histtype='stepfilled', alpha=0.3, bins=range(1,15))
        plt.hist(n_iso[1], label='IsoSeq', **kwargs)
        plt.hist(n_refseq_iso[1], label='RefSeq mRNA',**kwargs)
        plt.xlabel('Number of Isoforms per gene')
        plt.legend(loc='upper right')
        plt.yscale=('log')
        plt.savefig(out_fn)
        plt.close()
        kwargs['bins']=range(1,20)
        plt.hist(df['Splice_IsoSeq']/2+1, label='IsoSeq',normed=True, **kwargs)
        plt.hist(refseq_df['exonCount'], label='RefSeq mRNA',normed=True,**kwargs)
        out_fn=out_fn=os.path.splitext(fn)[0]+'_nExons.png'
        plt.xlabel('Number of exons per gene')
        plt.legend(loc='upper right')
        #plt.yscale=('log')
        plt.savefig(out_fn)
        plt.close()

if __name__=='__main__':
    css_stats, css_stats_group=get_stats(dir, "_ccs_report.txt", sumGroups=True)
    css_stats.iloc[0]/1e6
    css_stats_group.iloc[1]/css_stats_group.iloc[0]
    css_stats_group['Treatment']+css_stats_group['Control']


    lima_stats,lima_stats_group=get_stats(dir, "_demux.lima.summary", sumGroups=True)
    lima_stats_group.iloc[1]/lima_stats_group.iloc[0]
    lima_stats_group.loc['Undesired 3p--3p pairs']/lima_stats_group.iloc[0]
    lima_stats_group['Treatment']+lima_stats_group['Control']
