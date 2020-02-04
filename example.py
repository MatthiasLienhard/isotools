import isotools
import pandas as pd
from  isotools.transcriptome import Transcriptome
from importlib import reload
import copy
from statistics import mean
import logging
#reload(isotools.transcriptome)
import pickle
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns
import pomegranate #package to fit mixtures
import numpy as np
out_path='/project/42/pacbio/hecatos_isoseq'

isoseq_bam_fn='/project/42/pacbio/hecatos_isoseq/05-align/all_aligned_200129.sorted.bam'
isoseq_bam_fn='/project/42/pacbio/hecatos_isoseq/05-align/Control_aligned.sorted.bam'

#isoseq_fn='/project/42/pacbio/hecatos_isoseq/06-collapse/all_isoseq_collapsed.gff'
refseq_fn='/project/42/references/refseq/RefSeq_GRCh38_20181116_sorted.gff.gz'
#refseq_fn='refseq_head.gff.gz'
ensembl_fn='/project/42/references/Ensemblv98/Homo_sapiens.GRCh38.98.sorted.gff3.gz'
genome_fn='/project/42/references/GRCh38/hg38_chr_only.fa'
out_fn=isoseq_bam_fn.replace('05-align', '06-isotools').replace('.sorted.bam', '_isotools')
#reference=isotools.transcriptome.import_transcripts(ensembl_fn)
#with open('ensembl.pkl', 'wb') as f:
#    pickle.dump(reference, f, pickle.HIGHEST_PROTOCOL)
#with open('ensembl.pkl', 'rb') as f:
#    reference=pickle.load(f)
with open('refseq.pkl', 'rb') as f:
    reference=pickle.load(f)
#todo: remove/mark fusion genes in reference

chrom=[str(i+1)for i in range(22)]+['X','Y']
groups=[['m54070_190315_132401','m54070_190318_123927','m54070_190319_090034','m54070_190320_052230','m64080_200120_120529'], #
        ['m54070_190316_094451','m54070_190321_014433','m54070_190321_220654','m54070_190322_182901','m64080_200120_120529']]  #treatmentm64080_200121_181850
isoseq=Transcriptome(pacbio_fn=isoseq_bam_fn,ref=reference, chromosomes=chrom, groups=groups)
isoseq.find_truncations()
isoseq.find_biases(genome_fn) #this looks up junction types, direct repeats at junctions and downstream genomic a content

isoseq.add_filters()

sum(g.n_transcripts for g in isoseq)#367321
sum(g.n_transcripts-(len(g.data['truncated5']) if 'truncated5' in  g.data else 0) for g in isoseq)#276633

th=3
n_tr=[sum(1 for g in isoseq for tr in g.transcripts if tr['nZMW']>=th) for th in range(2,20)]
#[108065, 60834, 42257, 32714, 26986, 23193, 20379, 18250, 16592, 15230, 14137, 13183, 12395, 11723, 11114, 10544, 10075, 9665]
#[367321, 218758, 156876, 123454, 102438, 88023, 77359, 69411, 62962, 57817, 53494, 49901, 46882, 44270, 41916, 39897, 38136, 36341]
g=next(iter(isoseq))
g.to_gtf()
date='20200203'
df=isoseq.transcript_table(extra_columns=['length','n_exons','exon_starts','exon_ends' ,'all_canonical', 'grouped_nZMW','direct_repeat_len','template_switching','downstream_A_content','alt_splice','junction_type','truncation','filter'])
df.to_csv(f'{out_path}/tables/isoseq_transcripts_{date}.table', sep='\t',quoting=3, index=False) #3==QUOTE_NONE
g_ids={g.id:[g.name] for g in isoseq}
isoseq.write_gtf(f'{out_path}/tables/isoseq_transcripts_filtered_{date}_test.gtf', use_gene_name=True, remove={'a_th':.5, 'truncation':True, 'rrts':True})
isoseq.write_gtf(f'{out_path}/tables/isoseq_transcripts_filtered_novel_{date}.gtf', use_gene_name=True, include={'novel':True}, remove={'a_th':.5, 'truncation':True, 'rrts':True})

trunc=(g for g in isoseq if 'truncated5' in  g.data)

g=next(trunc)
print(g)
print(g.data['truncated5'])
for i,tr in enumerate(g.transcripts):
   print('{} - {}'.format(i,tr['exons']))
for read in align:

#chimeric reads
from pysam import  AlignmentFile
import re 
from isotools.transcriptome import junctions_from_cigar, overlap

cigarid={c:i for i,c in enumerate('MIDNSHP=XB')}
align=AlignmentFile(isoseq.pacbio_fn, 'rb')
for read in align:
    tags=dict(read.tags)
    if read.reference_name == 'MT' or tags['is']<10:
        continue
    if 'SA' in tags and not read.flag & 0x800:
        strand='-' if read.is_reverse else '-'
        exons_1=junctions_from_cigar(read.cigartuples, read.reference_start)
        snd_cigar=tags['SA'].split(',')
        snd_cigartuplestring=re.findall(r'(\d+)(\w+?)', snd_cigar[3])
        snd_cigartuple=[(cigarid[c],int(n)) for n,c in snd_cigartuplestring]
        exons_2=junctions_from_cigar(snd_cigartuple, int(snd_cigar[1]))
        dist=min(abs(exons_1[0][0]-exons_2[-1][1]),abs(exons_2[0][0]-exons_1[-1][1]))
        if read.reference_name == snd_cigar[0] and overlap((exons_1[0][0],exons_1[-1][1]),(exons_2[0][0],exons_2[-1][1])):
            continue
        print(f"{read.qname}, len={len(read.seq)}\t({tags['is']})\tchr{read.reference_name}:{read.reference_start}-{read.reference_end}({strand}) len={sum([e[1]-e[0] for e in exons_1])} [{dist}]\tchr{snd_cigar[0]}:{snd_cigar[1]}-{exons_2[-1][1]}({snd_cigar[2]}) len={sum([e[1]-e[0] for e in exons_2])}")
        #break


read.qname

roi='transcript/90154'
reads=[]
align=AlignmentFile(isoseq.pacbio_fn, 'rb')
for read in align:
    if read.qname == roi:
        reads.append(read)

for r in reads:
    print(r)

#plot A content

acontent=np.array(pd.to_numeric(df['downstream_A_content'])).reshape(-1,1)*30
model = pomegranate.gmm.GeneralMixtureModel.from_samples(pomegranate.distributions.PoissonDistribution, 2, acontent)
x = np.arange(0, max(acontent)+.01,0.05)
plt.figure(figsize=(10, 5))
plt.hist(acontent, bins=29, normed=True)
y0=model.distributions[0].probability(x)*np.exp(model.weights[0])
y1=model.distributions[1].probability(x)*np.exp(model.weights[1])
i_idx=np.argwhere(np.diff(np.sign(y0 - y1))).flatten() #index of intersection
plt.plot(x, y0, label="lambda={:.3f}, w={:.3f}".format(model.distributions[0].parameters[0],np.exp(model.weights[0])))
plt.plot(x, y1, label="lambda={:.3f}, w={:.3f}".format(model.distributions[1].parameters[0],np.exp(model.weights[1])))
plt.plot(x, model.probability(x), label="Mixture")
plt.axvline(x=x[i_idx][0], label='intersection:{:.1f}'.format(x[i_idx][0]), c='k')
plt.legend(fontsize=14)
plt.title("downstream A content", fontsize=14)
plt.ylabel("Probability Density", fontsize=14)
plt.savefig(plot_path+'/downstream_a_content.png')
# how many are above?
sum(pd.to_numeric(df['downstream_A_content'])<.5)/len(df) #18% of transcripts have more than 50% A
ne_normal=df.loc[pd.to_numeric(df['downstream_A_content'])<.5,'n_exons']
ne_highA=df.loc[pd.to_numeric(df['downstream_A_content'])>.5,'n_exons']
sns.distplot( ne_normal , color="skyblue", label="other")
sns.distplot( ne_highA , color="red", label="high genomic A downstream")
plt.legend()
plt.show() #twice as many have one exon only

#cases that might be template switching, since both sj are not in catalog
sel=list()
for g in isoseq:
    for tnr,tr in enumerate(g.transcripts):
        if tr['support'] is not None:
            for a in tr['support']['novel_sj2']:
                if a+1 in tr['support']['novel_sj1']:
                    sel.append((a, tnr,g))

sel=list()
for g in isoseq:
    sg=isotools.splice_graph.SpliceGraph(g)
    candidates=list(sg.find_ts_candidates())
    for st, end, js, ls, tr_idx in candidates:
        if ls>0 and ls/js>2:
            print('candidate {}:{}-{} {}/{}'.format(g.chrom,st-10, end+10, js, ls))
            sel.append(g)



len(sel) #6418 transcripts
from Bio import pairwise2

for j, (i, trn,g) in enumerate(sel):
    tr=g.transcripts[trn]
    print('{}:{}-{} {}'.format(tr['chr'],tr['exons'][i][1]-10,tr['exons'][i+1][0]+10, tr['junction_type'][i] ))
    print('ownscore={}'.format(tr['ts_score'][i][0]))
    seq=tr['ts_score'][i][1]
    sg=isotools.splice_graph.SpliceGraph(g)
    candidates=list(sg.find_ts_candidates())
    for  st, end, js, ls, tr_idx in candidates:

        print('candidate {}:{}-{} {}/{}'.format(g.chrom,st-10, end+10, js, ls))
    #print( pairwise2.format_alignment(*pairwise2.align.globalms(*seq,1,-1,-10,-1)[0]   ) )
    if j>10:
        break

import isotools.splice_graph
g=sel[5][2]
sg=isotools.splice_graph.SpliceGraph(g)
weights=[tr['nZMW'] for tr in g.transcripts]
#[tr['ts_score'][i-1]



#junction type statistics
jt=[t for g in isoseq for tr in g.transcripts for t in tr['junction_type']]
jt_known=[t for g in isoseq for tr in g.transcripts for i,t in enumerate(tr['junction_type']) if tr['support'] is not None and i-1 not in tr['support']['novel_sj1'] and i not in tr['support']['novel_sj2'] ]
jt_new=[t for g in isoseq for tr in g.transcripts for i,t in enumerate(tr['junction_type']) if tr['support'] is None or i-1 in tr['support']['novel_sj1'] or i in tr['support']['novel_sj2'] ]
jtc=Counter(jt) #97.77% GT-AG
jtc_known=Counter(jt_known)
jtc_new=Counter(jt_new)

for t,c in jtc.most_common(20):
    print('all {}: {:.4f}%'.format(t,c*100/len(jt)))
for t,c in jtc_known.most_common(20):
    print('known {}: {:.4f}%'.format(t,c*100/len(jt_known)))
for t,c in jtc_new.most_common(20):
    print('new {}: {:.4f}%'.format(t,c*100/len(jt_new)))

#get a noncanonical gene
for g in isoseq:
    if any(t != 'GTAG'  for tr in g.transcripts for t in tr['junction_type']):
        break





#get a gene with "gapped exons"
for g in isoseq:
    if any('gapped_exon' in tr['support']['sType'] for tr in g.transcripts if tr['support'] is not None):
        break

#{'strand': '-', 'chr': '17', 'ID': 'transcript/67090', 'exons': [[49404053, 49404977], [49405028, 49405204], [49406758, 49406854], [49409041, 49409158], [49409330, 49409473], [49411677, 49411839], [49413176, 49413291], [49414834, 49414871]], 'source': ['transcript/67090'], 'nZMW': 2, 'support': {'ref_gene_name': 'PHB', 'ref_gene_id': 'gene:ENSG00000167085', 'ref_transcript': 'transcript:ENST00000419140', 'ref_tss': 49414884, 'ref_pas': 49405072, 'ref_len': 907, 'ref_nSJ': 14, 'exI': 802, 'sjI': 12, 'exIoU': 0.42773333333333335, 'sjIoU': 0.75, 'sType': {'alternative_polyA': ['49404053-49404977'], 'exon_skipping': ['49413291~49414834']}}, 'ts_score': [6, 3, 1, 2, 3, 1, 5], 'junction_type': ['AGCC', 'GTAG', 'GTAG', 'GTAG', 'GTAG', 'GTAG', 'GTAG'], 'downstream_A_content': 0.1}
exons= [[49404053, 49404977], [49405028, 49405204], [49406758, 49406854], [49409041, 49409158], [49409330, 49409473], [49411677, 49411839], [49413176, 49413291], [49414834, 49414871]]
def get_gene(reference, id):
    for chrom,t in reference.items():
        for g in t:
            if g.name==id or g.id==id:
                return(g)
            for tr in g.transcripts:
                if tr['ID']==id:
                    return(tr)

tr=get_gene(reference, 'transcript:ENST00000419140')



gi=iter(isoseq)
g=next(gi)

ts_score_known=np.array([s for g in isoseq for tr in g.transcripts if tr['support'] is not None for s in tr['ts_score']])
ts_score_unknown=np.array([s for g in isoseq for tr in g.transcripts if tr['support'] is None for s in tr['ts_score']])
plt.figure(figsize=(10, 5))
plt.hist(ts_score_unknown, bins=20, normed=True)
plt.hist(ts_score_known, bins=20, normed=True)

plt.show()

list(g.transcripts[0].keys())
g.transcripts[0]['ts_score'] #this requires refinement
g.transcripts[0]['junction_type'] #GTAG is canonical
g.transcripts[0]['downstream_A_content'] # what threshold is alarming? 50%
#whats the optimal length? look at the distribution and some examples to see how much is required (currently default 30 bases)





g1=get_gene('gene28069')
g2=get_gene('gene28070')
new=[[47541308, 47542306], [47543541, 47543663], [47544388, 47544529], [47545431, 47545620], [47549304, 47549384], [47550657, 47550847], [47551743, 47551952], [47553478, 47553506]]
novel={'10':[{'chr':10, 'ID':'test', 'exons':new, 'strand':'-'}]}
isoseq.add_novel_transcripts(novel)

chrom='10'
tr={'chr':10, 'ID':'test', 'exons':new, 'strand':'-'}
candidates=isoseq.data[chrom].overlap(tr['exons'][0][0], tr['exons'][-1][1])
merge_genes=set()
for ol_gene in candidates:
    if ol_gene.strand != tr['strand']:
        continue
    if any(isotools.transcriptome.is_same_gene(tr['exons'], oltr['exons'],spj_iou_th, reg_iou_th) for oltr in ol_gene.transcripts):
        merge_genes.add(ol_gene)
        print(ol_gene)

g=next(iter(g for g in isoseq if toi in t['source']  ))
g.transcripts[toi]

print(isoseq)
isoseq.collapse_transcripts(fuzzy_junction=0) 
print(isoseq)
#isoseq.collapse_transcripts(fuzzy_junction=5) 
#print(isoseq)
isoseq_bare=copy.deepcopy(isoseq)
trnames_bare={tr for g in isoseq_bare for tr in g.data['transcripts']}
#isoseq.write_gtf(isoseq_bam_fn.replace('.bam','.gtf'))
refseq=Transcriptome(refseq_fn, chromosomes=isoseq.chromosomes)

#print(refseq)
isoseq.collapse_to_reference(refseq)
print(isoseq)

#ensembl=Transcriptome(ensembl_fn, chromosomes=isoseq.chromosomes)

#pb8tr=next(iter(isoseq['PB_8'])).transcripts
#genes={g for g in isoseq_bare if next(iter(g.transcripts.keys())) in pb8tr}
#isotools.utils.get_support(next(iter(next(iter(genes)).transcripts.values()))['exons'], refseq.data, '1', False)
isoseq.write_gtf(out_fn+'_annotated.gtf', source='pacbio')


#refseq_collapsed=copy.deepcopy(refseq)
#refseq_collapsed.collapse_transcripts()


gtab=isoseq.gene_table().sort_values(['chr', 'begin', 'end'])
#print(isoseq_collapsed)
gtab.to_csv(out_fn+'_genes.txt', sep='\t',index=False)
#ttab=isoseq.transcript_table(region='1:9000000-10000000', extra_columns=['length','alt_splice'])

ttab=isoseq.transcript_table( extra_columns=['length','n_exons','alt_splice']).sort_values(['chr', 'begin', 'end'])
ttab.to_csv(out_fn+'_transcripts.txt', sep='\t',index=False)


import stats
stats.plot_altsplice(ttab,out_fn=out_fn+'_altsplice.png' )


# altsplice stats:
from pysam import TabixFile, AlignmentFile, FastaFile
from Bio import pairwise2
genome_fn='/project/42/references/GRCh38/hg38_chr_only.fa'
genome_fh=FastaFile(genome_fn)


g=next(iter(isoseq))
g.find_threeprime_a_content(genome)
tr=g.transcripts[0]
delta=10
introns=[(tr['exons'][i][1], tr['exons'][i+1][0]) for i in range(len(tr['exons'])-1)]
intron_seq=[[genome.fetch(g.chrom, pos-delta, pos+delta) for pos in i] for i in introns]
align=[[a==b for a,b in zip(*seq)] for seq in intron_seq ]





typeOI='gapped_exon'
gaps=list()
#look at direct repeats and template switching
for g in (g for g in isoseq if any(typeOI in t['support']['sType'] for t in g.transcripts.values() if t['support'] is not None)):
    for tid,t in ((tid, t) for tid, t in  g.transcripts.items() if typeOI in t['support']['sType']):
        print(t['support']['sType'][typeOI])
        for gap in t['support']['sType'][typeOI][0].split('-')[1:-1]:
            intron=[int(i) for i in gap.split('~')]
            print(intron[1]-intron[0])            
            #rl=find_ts_repeat_length(seq)
            gaps.append((intron[1]-intron[0],check_rt_switching(g.chr , intron, genome, min_score=8 )))
            #align=pairwise2.align.localms(seq[0], seq[1], 1,-2, -2, -1)[0]
            #print(pairwise2.format_alignment(*align))            
            #print(align)
            #gaps.append((intron[1]-intron[0],align[2] ))

sum(g[1] for g in gaps)/len(gaps)
#12,8% template switching

ctl=list()
for i,g in enumerate(refseq):
    if i>300:
        break
    for tid,t in g.transcripts.items():
        for jid in range(len(t['exons'])-1):
            intron=(t['exons'][jid][1],t['exons'][jid+1][0])
            ctl.append((intron[1]-intron[0],check_rt_switching(g.chr , intron, genome, min_score=8 )))


sum(g[1] for g in ctl)/len(ctl)
#0.4%




def check_rt_switching(chr, pos, genome,wiggle=1, min_score=4):
    delta=wiggle+min_score
    seq=[genome.fetch(g.chr, i-delta, i+delta) for i in pos]
    align=pairwise2.align.localms(seq[0], seq[1], 1,-2, -2, -1)
    #wig=[a[0][:a[3]].count('-')-a[1][:a[3]].count('-') for a in align]
    try:
        used_align=next(iter(a for a in align if a[2]>=min_score and abs(a[0][:a[3]].count('-')-a[1][:a[3]].count('-')) <= wiggle))
        print(pairwise2.format_alignment(*used_align))    
        return True
    except StopIteration:
        return False



 if True:   
    cmp=[c1==c2 for c1,c2 in zip(*seq)]
    sum=0
    found=False
    for i,eq in enumerate(cmp):
        if i==10:
            found=True
        if eq:
            sum+=1
        else:
            if not found:
                sum=0
            else:
                return(sum)
    return 0


#look for toxic /interesting cases
typeOI='gapped_exon'
g=next(iter(g for g in isoseq if any(typeOI in t['support']['sType'] for t in g.transcripts.values() if t['support'] is not None)))
alt_t=next(iter(t for t in  g.transcripts.values() if typeOI in t['support']['sType']))
ref_t=next(iter(refseq[alt_t['support']['ref_gene']])).transcripts[alt_t['support']['ref_transcript']]
isotools.utils.get_splice_type(ref_t['exons'], alt_t['exons'], g.strand=='-')




#some stats:
sum(ttab['n_exons']>1)
#50.910 multi exon transcripts

sum(ttab['splice_type']=='NA')
len(ttab['gene_name'].loc[ttab['splice_type']=='NA'].unique())
sum(ttab['n_exons'].loc[ttab['splice_type']=='NA']>1)
len(ttab['gene_name'].loc[(ttab['n_exons']>1) & (ttab['splice_type']=='NA')].unique())
#5.636 transcripts from 5.240 genes are novel, 1.465 > 1 exon

identical=[True if 'splice_identical' in alt else False for alt in  ttab['splice_type']]
sum(identical)
len(set(ttab['gene_name'].loc[identical]))
#7.710 isoforms “splice identical” (splice IoU=1) to Refseq (5.510 genes)

#gapped exons
gapped_exons=[v for t in ttab['splice_type'] if isinstance(t,dict) for k,v in t.items() if k == 'gapped_exon']
for g in isoseq:
    for tid,t in g.transcripts.items():
        if t['support'] is not None:
            if 'gapped_exon' in t['support']['sType']:
                print ('{} {}: {}'.format(tid,t['exons'],t['support']['sType']['gapped_exon'] ))
                raise Exception
            #if 'gapped_exon' in t.support

#SJIOU>1?
strange=[v>1 if v != 'NA' else False for v in  ttab['splice_IoU']]
#not defined splice type
ttab.loc[ttab['splice_type']=={}]
trids=set(ttab.loc[ttab['splice_type']=={}]['transcript_name'])
tr=[v for g in isoseq for t,v in g.transcripts.items() if t in trids]
[len(t['exons']) for t in tr]
#all are length 1

[t for g in [g for g in isoseq if any(t in trids for t in g.transcripts)] for t in g if t in trids]
[t for g in isoseq for t in g.transcripts if t in trids]
toi='transcript/108245'

for toi in trids:
    g=next(iter(g for g in isoseq if toi in g.transcripts))
    #g.transcripts[toi]['support']
    #support=isotools.utils.get_support(g.transcripts[toi]['exons'],refseq.data, g.chr, g.strand=='-'  )
    alt_exons=g.transcripts[toi]['exons']
    refseq_oi=g.transcripts[toi]['support']['ref_transcript']
    ref_exons=[v for g in refseq for t,v in g.transcripts.items() if t ==refseq_oi][0]['exons']
    print('{}: {}'.format(g,isotools.utils.get_relation(ref_exons, alt_exons))



refseq_oi='NM_018085.4'
next(iter(g for g in refseq if refseq_oi in g.transcripts)).transcripts[refseq_oi]['exons']

print(isoseq)
print(next(iter(isoseq)))

isoseq.history

print(isoseq)



def get_splice_type(ref, alt, is_reversed=False):
    if len(ref) == 0:
        return(['novel/unknown'])
    if len(alt) ==1:
        return['unspliced_fragment']
    types = ['alternative_donor', 'alternative_acceptor', 'alternative_promoter', 'alternative_polyA',
             'truncated5', 'truncated3', 'exon_skipping', 'novel_exon', 'gapped_exon', 'retained_intron']
    types = {t: [] for t in types}
    relation = isotools.utils.get_relation(ref, alt)
    # alt exons that do not overlap any ref
    first = next(i for i, rel in enumerate(relation) if rel) #first ref exon that overlaps an alt exon
    last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
    novel = set(range(len(alt)))-{r[0] for x in relation for r in x}
    if 0 in novel:
        types['alternative_polyA' if is_reversed else 'alternative_promoter'].append(
            '{}-{}'.format(*alt[0]))
    if len(alt)-1 in novel:
        types['alternative_promoter' if is_reversed else 'alternative_polyA'].append(
            '{}-{}'.format(*alt[-1]))
    if len(novel - {0, len(alt)-1}) > 0:
        # todo: distinguish completely novel from novel for this ref_transcript
        types['novel_exon'] += ['{}-{}'.format(*e)
                                for i, e in enumerate(alt) if i in novel-{0, len(alt)-1}]
    # if 0 not in novel and not types['novel_exon'] and len(relation[0]) == 0:
    all_match = (len(novel) == 0 and #no novel alt exons
        all([len(rel) == 1 for rel in relation[first:(last+1)]]) and #each alt exon corresponds to exactly one ref exon
        all([rel[0][1] == 3 for rel in relation[(first+1):last]]) and #all but the first and last alt exons are the same
        (relation[first][0][1] & 1) and #the first exon has matching second splice site
        (relation[last][0][1] & 2)) #the last exon has matching first splice site
    if first > 0 and all_match:
        types['truncated3' if is_reversed else 'truncated5'].append(
            '{}-{}'.format(*alt[0]))
    # if len(alt)-1 not in novel and not types['novel_exon'] and len(relation[-1]) == 0:
    if last < len(relation)-1 and all_match:
        types['truncated5' if is_reversed else 'truncated3'].append(
            '{}-{}'.format(*alt[-1]))
    for i in range(first, last+1):
        rel = relation[i]
        if len(rel) > 1:  # more than one alt exon for a ref exon
            types['gapped_exon'].append(
                "~".join(['{}-{}'.format(*alt[alt_i]) for alt_i, splice in rel]))
        if rel and i > first and relation[i-1] and rel[0][0] == relation[i-1][-1][0]:
            types['retained_intron'].append('{}-{}'.format(*alt[rel[0][0]]))
        # fst splice site does not correspond to any alt exon
        if rel and i < last and rel[-1][1] & 1 == 0 and i < last:
            delta = alt[rel[-1][0]][1]-ref[i][1]
            types['alternative_donor' if is_reversed else 'alternative_acceptor'].append(
                (alt[rel[-1][0]][1], delta))
        if rel and rel[0][1] & 2 == 0 and i > first:
            delta = alt[rel[0][0]][0]-ref[i][0]
            types['alternative_acceptor' if is_reversed else 'alternative_donor'].append(
                (alt[rel[0][0]][0], delta))
        if not rel and i > first and i < last:  # exon skipping
            pre=next(((j,relation[j][-1]) for j in reversed(range(i)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            suc=next(((j,relation[j][0]) for j in range(i+1, len(relation)) if relation[j]),[]) #first ref exon that overlaps an alt exon
            # only predecessing as well as successing splice site is identical
            #if pre and suc and pre[1][1] & 1 and suc[1][1] & 2:
            types['exon_skipping'].append(
                    '{}~{}'.format(alt[pre[1][0]][1], alt[suc[1][0]][0]))
    # return only types that are present
    return {k: v for k, v in types.items() if v}
