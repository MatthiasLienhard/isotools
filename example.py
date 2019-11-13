import isotools
import pandas as pd
from  isotools.transcriptome import Transcriptome
from importlib import reload
from statistics import mean
#reload(isotools.transcriptome)
isoseq_bam_fn='/project/42/pacbio/hecatos_isoseq/05-minimap/all_merged_minimap.bam'
isoseq_fn='/project/42/pacbio/hecatos_isoseq/05-minimap/all_merged_collapsed.gtf.gz'
refseq_fn='/project/42/references/refseq/RefSeq_GRCh38_20181116_sorted.gff.gz'

isoseq_collapsed=Transcriptome(isoseq_fn,type='gtf', rescue_genes=True)
isoseq=Transcriptome(isoseq_bam_fn)
refseq=Transcriptome(refseq_fn)
isoseq.collapse_transcripts()
isoseq.add_splice_graphs(isoseq)
isoseq.add_reference_support(refseq)







#look for toxic /interesting cases
#genes with multiple isoforms equal to the same refseq gene (which should not be there)
toi=('splice_identical','truncated5', 'truncated3')
for chrom, tree in isoseq.items():
    for i,gene in enumerate(tree):
        tr_identical={trid:tr['support']['ref_transcript'] for trid,tr in gene.data['transcripts'].items() if any(term in tr['support']['sType'] for term in toi)}
        n_found={refseq:list(tr_identical.values()).count(refseq) for refseq in set(tr_identical.values())}
        try:
            multi=next(refseq for refseq,n in n_found.items() if n>1)
        except StopIteration:
            pass
        else: #there is more than one
            exons=[tr['exons'] for tr in gene.data['transcripts'].values() if any(term in tr['support']['sType'] for term in toi) and tr['support']['ref_transcript']==multi]
            raise Exception('chrom {} gene {}: {}'.format(chrom, i,multi))

#Exception: chrom 1 gene 106: NM_001351365.1
#list(isoseq['1'])[106]


for e in exons:
 print(e)

ref_genes=refseq[chrom].overlap(*gene[:2])
pd.DataFrame(isotools.splice_graph.compute_support(ref_genes, exons[0]))
ref_exons=next(g.data['transcripts'][multi]['exons'] for g in ref_genes if multi in g.data['transcripts'])

[isotools.splice_graph.get_splice_type(ref_exons, e, gene.data['strand']=='-') for e in exons]

relation=isotools.splice_graph.get_relation(ref_exons,exons[1])
isotools.splice_graph.get_intersects(ref_exons,exons[1])
isotools.splice_graph.is_truncation(ref_exons,exons[1], debug=True)
isotools.splice_graph.check_truncation(exons[0],exons[1])
isotools.splice_graph.is_truncation(exons[0], exons[1])
relation=isotools.splice_graph.get_relation(exons[0], exons[1])
any(len(r) > 1 for r in relation)#one exon from 1 has several corresponding exons from 2
not any(r for r in relation ) #no relation
first=next(iter([i for i,r in enumerate(relation) if r]))
last=len(relation)-next(iter([i for i,r in enumerate(reversed(relation)) if r]))-1
sum(1 for r in relation if r) != last-first+1  # some exons in 1 do not have corresponding exons in 2
first>0 and relation[first][0][0]>0
last < len(exons[0])-1 and relation[last][0][0]<len(exons[1])-1 
last!=first and ( #more than one exon
            (relation[first][0][1] & 1 ==0) or # 2nd splice junction must match
            (relation[last][0][1] & 2 == 0)) # fst splice junction must match
relation[first][0][0]!=first and (first==0 ) == (exons[0][first][0] >= exons[1][relation[first][0][0]][0]) # first exon of trunkated version is actually longer
relation[last][0][0] != last-first and (last==len(exons[0])) == (exons[0][last][1] < exons[1][relation[last][0][0]][1]) # last exon of trunkated version is actually longer
last-first > 3 and any(relation[i][0][1]!=3 for i in range(first+1, last)) #intermediate exons do not fit



collapsed=dict()
for tr1_name, tr1 in gene.data['transcripts'].items():
    for tr2_name, tr2 in collapsed.items():
        truncation=isotools.splice_graph.check_truncation(tr1['exons'], tr2['exons'])
        if truncation:
            #todo: check thresholds
            collapsed[tr2_name]=isotools.splice_graph.merge_transcripts(tr1,tr2)
            
            break
    else:
        collapsed[tr1_name]=tr1


gene=Interval(39738871, 39763913, {'strand': '+', 'Name': 'PB_813', 'source': {'PB_813_1': ['all_merged_transcript/78754|full_length_coverage=2|length=1724|num_subreads=41']}, 'transcripts': {'PB_813_1': {'exons': [[39738871, 39738923], [39740164, 39740260], [39741378, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753053], [39753285, 39754042]], 'source': ['all_merged_transcript/78754|full_length_coverage=2|length=1724|num_subreads=41']}, 'PB_813_2': {'exons': [[39738872, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39753705]], 'source': ['all_merged_transcript/91637|full_length_coverage=211|length=1287|num_subreads=60']}, 'PB_813_3': {'exons': [[39738878, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745413, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39753706]], 'source': ['all_merged_transcript/91851|full_length_coverage=14|length=1242|num_subreads=60']}, 'PB_813_4': {'exons': [[39738878, 39738923], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39763650, 39763913]], 'source': ['all_merged_transcript/96465|full_length_coverage=3|length=1120|num_subreads=60']}, 'PB_813_5': {'exons': [[39738881, 39738923], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39756762]], 'source': ['all_merged_transcript/11115|full_length_coverage=3|length=4328|num_subreads=60']}, 'PB_813_6': {'exons': [[39738881, 39738923], [39740164, 39740263], [39741365, 39741409], [39741894, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39753699]], 'source': ['all_merged_transcript/49163|full_length_coverage=6|length=2543|num_subreads=60']}, 'PB_813_7': {'exons': [[39738881, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39763688, 39763913]], 'source': ['all_merged_transcript/98168|full_length_coverage=26|length=1082|num_subreads=60']}, 'PB_813_8': {'exons': [[39738883, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39754042]], 'source': ['all_merged_transcript/81997|full_length_coverage=5|length=1613|num_subreads=60']}, 'PB_813_9': {'exons': [[39738883, 39738933], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753285, 39753698]], 'source': ['all_merged_transcript/89788|full_length_coverage=2|length=1312|num_subreads=17']}, 'PB_813_10': {'exons': [[39738883, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39762516, 39762651], [39763688, 39763913]], 'source': ['all_merged_transcript/95337|full_length_coverage=6|length=1218|num_subreads=60']}, 'PB_813_11': {'exons': [[39738891, 39738923], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39748902, 39749088], [39752909, 39753052], [39753286, 39753696]], 'source': ['all_merged_transcript/97763|full_length_coverage=2|length=1127|num_subreads=38']}, 'PB_813_12': {'exons': [[39738892, 39738931], [39740164, 39740263], [39741365, 39741409], [39741894, 39741921], [39743215, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39763688, 39763913]], 'source': ['all_merged_transcript/97510|full_length_coverage=2|length=1082|num_subreads=22']}, 'PB_813_13': {'exons': [[39743214, 39743297], [39743823, 39743924], [39745374, 39745498], [39748902, 39749088], [39752909, 39753052], [39753286, 39753706]], 'source': ['all_merged_transcript/98779|full_length_coverage=3|length=1061|num_subreads=60']}}})



exons1=[(i*20+10, (i+1)*20) for i in range(7)]
exons1=[(10, 20), (30, 40), (50, 60), (70, 80), (90, 95)]
exons2=[(35,40), (50, 60), (70, 80), (90, 100), (110, 120), (130, 140), (150, 160), (170, 180), (190, 200), (210, 220)]

isotools.splice_graph.check_truncation(exons1,exons2)


