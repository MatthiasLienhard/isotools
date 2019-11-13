import os
import warnings
import copy
from pysam import TabixFile, AlignmentFile
#from Bio import SeqIO, Seq, SeqRecord
from intervaltree import IntervalTree, Interval
from tqdm import tqdm
import isotools.utils

def import_transcripts(fn, type='auto', **kwargs):
    if type=='auto':        
        type=os.path.splitext(fn)[1].lstrip('.')
        if type=='gz':
            type=os.path.splitext(fn[:-3])[1].lstrip('.')
    if type == 'bam':
        genes= import_bam_transcripts(fn,  **kwargs)
    elif type == 'gtf':
        genes= import_gtf_transcripts(fn,  **kwargs)
    elif type == 'gff':
        genes= import_gff_transcripts(fn,  **kwargs)
    else:
        raise ValueError('unsupportet file type: {}'.format(type))
    
    return genes

def import_gtf_transcripts(fn, genes=None,chromosomes=None, rescue_genes=False):
    gtf = TabixFile(fn)   
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()
    if genes is None:
        genes=dict()   
    gids=set() 
    for line in tqdm(gtf.fetch(), smoothing=.1): 
        ls = line.split(sep="\t")
        
        if chromosomes is not None and ls[0] not in chromosomes:
            #warnings.warn('skipping line from chr '+ls[0])
            continue
        genes.setdefault(ls[0], IntervalTree())
        info = dict([pair.lstrip().split(' ', 1) for pair in ls[8].replace('"','').split(";") if pair])        
        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if ls[2] == "exon":
            # print(line)
            try:
                gtf_id = info['transcript_id']
                exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                warnings.warn("gtf format error: exon without transcript_id. Skipping line\n"+line)
        elif ls[2] == 'gene':
            if 'gene_id' not in info:
                warnings.warn("gtf format error: gene without gene_id. Skipping line\n"+line)
            else:
                info['strand'] = ls[6]
                info['chromosome'] = ls[0]
                genes[ls[0]][start:end] = info        
                gids.add(info['gene_id'])
        elif ls[2] == 'transcript':
            try:
                transcripts.setdefault(info["gene_id"], list()).append(info["transcript_id"])
            except KeyError:
                warnings.warn("gtf format error for transcript, skipping line\n"+line )
            else:
                if rescue_genes:
                    if info["gene_id"] in gids:
                        extend_gene(info["gene_id"],start, end, genes[ls[0]])
                    else:
                        info=copy.deepcopy(info)
                        del(info['transcript_id'])
                        info['strand'] = ls[6]
                        info['chromosome'] = ls[0]
                        genes[ls[0]][start:end] = info        
                        gids.add(info['gene_id'])                        
        else:
            skipped.add(ls[2])
    if len(skipped)>0:
        print('skipped the following categories: {}'.format(skipped))
    # sort the exons
    for tid in exons.keys():
        exons[tid].sort()
    # add transcripts to genes

    if not rescue_genes:
        #check for missing gene information
        missed_genes=0
        for gid in transcripts:
            if gid not in gids:
                missed_genes+=1
        if missed_genes:
            raise ValueError('no specific gene information for {}/{} genes'.format(missed_genes,missed_genes+sum((len(t) for t in genes.values()) )))
            


    for chrom in genes:
        for gene in genes[chrom]:
            g_id = gene.data['gene_id']
            t_ids = transcripts.get(g_id,  g_id)#if there are no transcripts use gene id
            for t_id in t_ids:
                try:
                    gene.data.setdefault('transcripts', dict())[t_id] = {'exons':exons[t_id]}
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    gene.data['transcripts'] = {t_id: {'exons':[tuple(gene[:2])]}}
                    
    return genes

def extend_gene(id,start, end, genes):
    for gene in genes.overlap(start, end):
        if gene.data['gene_id']==id:
            merged_gene=Interval(min(start,gene.begin), max(end, gene.end), gene.data)
            genes.add(merged_gene)
            genes.remove(gene)    
            return True
    return False


def import_gff_transcripts(fn, genes=None):
    # returns a dict interval trees for the genes, each containing the splice graph
    gff = TabixFile(fn)        
    chrom_ids = get_gff_chrom_dict(gff)
    exons = dict()  # transcript id -> exons
    transcripts = dict()  # gene_id -> transcripts
    skipped = set()    
    if genes is None:
        genes=dict()
    #takes quite some time... add a progress bar?
    for line in tqdm(gff.fetch(), smoothing=.1):  # parser=pysam.asGTF()): parser is strange and not helpfull
        ls = line.split(sep="\t")
        if ls[0] not in chrom_ids:
            #warnings.warn('unknown chromosome id :'+ls[0])
            # this happens at all NW/NT scaffolds
            continue
        chrom = chrom_ids[ls[0]]
        genes.setdefault(chrom, IntervalTree())
        try:
            info = dict([pair.split('=', 1) for pair in ls[8].split(";")])
        except ValueError:
            warnings.warn("GFF format error in infos (should be ; seperated key=value pairs). Skipping line: "+line)

        start, end = [int(i) for i in ls[3:5]]
        start -= 1  # to make 0 based
        if 'Name' not in info:
            for alt in 'standard_name', 'gene', 'product':
                try:
                    info['Name'] = info[alt]
                    break
                except KeyError:
                    pass
        if ls[2] == "exon":
            # print(line)
            try:
                gtf_id = info['Parent']
                exons.setdefault(gtf_id, list()).append((start, end))
            except KeyError:  # should not happen if GTF is OK
                warnings.warn("GFF format error: no Parent found for exon. Skipping line: "+line)
        # elif ls[2] in ['gene', 'pseudogene']:
        elif ls[2] == 'gene':
            info['strand'] = ls[6]
            info['chromosome'] = chrom
            genes[chrom][start:end] = info
        # those denote transcripts
        elif all([v in info for v in ['Parent', "ID", 'Name']]) and info['Parent'].startswith('gene'):
            transcripts.setdefault(info["Parent"], list()).append(
                (info["Name"], info["ID"]))
        else:
            skipped.add(ls[2])
    #print('skipped the following categories: {}'.format(skipped))
    # sort the exons
    for tid in exons.keys():
        exons[tid].sort()
    # add transcripts to genes
    for chrom in genes:
        for gene in genes[chrom]:
            g_id = gene.data['ID']
            t_ids = transcripts.get(g_id, [(gene.data['Name'], g_id)])
            for t_name, t_id in t_ids:
                try:
                    gene.data.setdefault('transcripts', dict())[t_name] = {'exons':exons[t_id]}
                except KeyError:
                    # genes without transcripts get a single exons transcript
                    gene.data['transcripts'] = {t_name: {'exons':[tuple(gene[:2])]}}
                    
    return genes

def import_bam_transcripts(bam_fn, genes=None, n_lines=None):
    gene_prefix='gene_'
    genes = dict()
    gene_number=0
    merged_number=0
    with AlignmentFile(bam_fn, "rb") as align:        
        stats = align.get_index_statistics()
        # try catch if sam/ no index
        total_reads = sum([s.mapped for s in stats])
        if n_lines is not None:
            total_reads=n_lines
        for i, read in tqdm(enumerate(align.fetch()), total=total_reads, unit='transcripts'):
            if read.flag & 0x800:
                continue #todo: include chimeric reads!
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            genes.setdefault(chrom, IntervalTree())
            exons = isotools.utils.junctions_from_read(read)
            pos = exons[0][0], exons[-1][1]
            found=None
            for gene in genes[chrom].overlap(*pos):                
                if gene.data['strand'] == strand and any([isotools.utils.is_same_gene(exons, tr['exons']) 
                            for tr in gene.data['transcripts'].values()]):
                    if found is None:
                        tr_name='{}_{}'.format(gene.data['Name'],len(gene.data['transcripts'])+1)
                        gene.data['transcripts'][tr_name]={'exons':exons, 'source':[read.query_name]}
                        if pos[0]<gene.begin or pos[1]>gene.end:
                            merged_gene=Interval(min(pos[0],gene.begin), max(pos[1], gene.end), gene.data)
                            genes[chrom].add(merged_gene)
                            genes[chrom].remove(gene)    
                            found=merged_gene
                        else:
                            found=gene
                    else:
                        #merge (current) gene and (already) found_gene, as both overlap with new transcript
                        #print('this happens!!')
                        merged_number+=1
                        found.data['transcripts'].update(gene.data['transcripts'])
                        if gene.begin<found.begin or gene.end > found.end:
                            merged_gene=Interval(min(found.begin,gene.begin), max(found.end, gene.end), found.data)
                            genes[chrom].add(merged_gene)
                            genes[chrom].remove(found)    
                            found=merged_gene
                        genes[chrom].remove(gene)                            
            if found is None: #new gene
                gene_number+=1
                gname=gene_prefix+str(gene_number)
                info={'strand':strand, 'Name':gname,'source':{gname+'_1':[read.query_name]}, 
                        'transcripts':{gname+'_1':{'exons':exons,'source':[read.query_name]}}}
                genes[chrom].addi(*pos, info)
            if n_lines==i:
                break
        
    return genes


def get_gff_chrom_dict(gff):
    # fetch chromosome ids
    chrom = {}
    for c in gff.contigs:
        #print ("---"+c)
        for line in gff.fetch(c, 1, 2):
            if line[1] == "C":
                ls = line.split(sep="\t")
                if ls[2] == "region":
                    info = dict([pair.split("=")
                                    for pair in ls[8].split(";")])
                    if "chromosome" in info.keys():
                        chrom[ls[0]] = info["chromosome"]
                        break
                        # print(line)
        else:
            chrom[c]=c # no specific regions entrie - no aliases
    return(chrom)



def collapse_isoforms(bam_fn='/project/42/pacbio/hecatos_isoseq/05-gmap/all_merged_hq_isoforms.bam'):
    genes = dict()
    gene_number = 0
    isoform_number = 0
    with AlignmentFile(bam_fn, "rb") as align:
        stats = align.get_index_statistics()
        # try catch if sam/ no index
        total_reads = sum([s.mapped for s in stats])
        for i, read in tqdm(enumerate(align.fetch()), total=total_reads):
            tr_name = read.query_name[:read.query_name.find('|')]
            chr = read.reference_name
            strand = '-' if read.is_reverse else '+'
            genes.setdefault(chr, IntervalTree())
            exons = isotools.utils.junctions_from_read(read)
            pos = exons[0][0], exons[-1][1]
            for gene in genes[chr].overlap(*pos):
                same_gene = False
                if gene.data[0] == strand:
                    for j, tr in enumerate(gene.data[1]):
                        if isotools.collapse.splice_identical(tr, exons): 
                            tr[0][0] = min(pos[0], tr[0][0])
                            tr[-1][1] = max(exons[-1][1], pos[1])
                            gene.data[3][j].append(pos)
                            same_gene = True
                            break
                        # same gene: at least one shared splice site
                        elif any([r[1] > 0 for x in isotools.utils.get_relation(exons, tr) for r in x if r]):
                            same_gene = True  # continue looking for a better match
                            
                    else:
                        if same_gene:
                            gene.data[1].append(exons)
                            gene.data[3].append([pos])
                            isoform_number += 1
                if same_gene:
                    break
            else:  # new gene
                genes[chr].addi(
                    *pos, (strand, [exons], 'PB_{}'.format(gene_number), [[pos]]))
                gene_number += 1
                isoform_number += 1



def get_read_sequence(bam_fn,reg=None, name=None):
    with AlignmentFile(bam_fn, "rb") as align:
        for i, read in enumerate(align.fetch(region=reg)):
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            exons = isotools.utils.junctions_from_read(read)
            if name is None or name == read.query_name:
                yield  read.query_name, read.query_sequence, chrom,strand, exons
