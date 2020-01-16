import isotools
import pandas as pd
from  isotools.transcriptome import Transcriptome
from importlib import reload
import copy
from statistics import mean
import logging
#reload(isotools.transcriptome)
import pickle



isoseq_bam_fn='/project/42/pacbio/hecatos_isoseq/05-align/all_aligned.sorted.bam'
#isoseq_fn='/project/42/pacbio/hecatos_isoseq/06-collapse/all_isoseq_collapsed.gff'
refseq_fn='/project/42/references/refseq/RefSeq_GRCh38_20181116_sorted.gff.gz'
#refseq_fn='refseq_head.gff.gz'
ensembl_fn='/project/42/references/Ensemblv98/Homo_sapiens.GRCh38.98.sorted.gff3.gz'
genome_fn='/project/42/references/GRCh38/hg38_chr_only.fa'
out_fn=isoseq_bam_fn.replace('05-align', '06-isotools').replace('.sorted.bam', '_isotools')
#reference=isotools.transcriptome.import_transcripts(refseq_fn)
#with open('refseq.pkl', 'wb') as f:
#    pickle.dump(reference, f, pickle.HIGHEST_PROTOCOL)

with open('refseq.pkl', 'rb') as f:
    reference=pickle.load(f)


isoseq=Transcriptome(pacbio_fn=isoseq_bam_fn,ref=reference)
print(isoseq)
isoseq.remove_chromosome('MT')
print(isoseq)
def get_gene(gid):
    for chrom, tree in reference.items():
        for g in tree:
            if g.id == gid:
                return(g)

g1=get_gene('gene28069')
g2=get_gene('gene28070')
new=[[47541308, 47542306], [47543541, 47543663], [47544388, 47544529], [47545431, 47545620], [47549304, 47549384], [47550657, 47550847], [47551743, 47551952], [47553478, 47553506]]

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
genome=FastaFile(genome_fn)
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
