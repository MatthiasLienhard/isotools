import logging
import gzip
import pandas as pd
from pysam import FastaFile
import pyhmmer
from tqdm import tqdm
import requests
import time
from ._utils import genomic_position, has_overlap


logger = logging.getLogger('isotools')


def parse_hmmer_metadata(fn):
    with gzip.open(fn) as f:
        entries = []
        entry = {}
        while (True):
            try:
                line = next(f).decode().strip()
            except StopIteration:
                metadata = pd.DataFrame(entries).set_index('AC')
                return metadata
            if line == '//':
                entries.append(entry.copy())
            elif line.startswith('#=GF'):
                _, k, v = line.split(maxsplit=2)
                entry[k] = v


def add_domains_to_table(table, transcriptome, source='annotation', categories=None, id_col='gene_id', modes=['trA-trB', 'trB-trA'],
                         naming='id', overlap_only=False, insert_after=None, **filter_kwargs):
    '''add domain annotation to table.

    :param table: A table, for which domains are derived.
        It should have at least one column with a gene id and one with a list of transcripts.
    :param transcriptome: The isotools.Transcriptome object, from which the domains are sourced.
    :param source: One of the three possible annotation sources.
    :param categories: A list of annotation categories, which are considered.
        If set to None (default) all annotation categories are considered.
    :param modes: A list of strings with column names of transcript ids used to look up domains.
        Multible columns can be combined to one mode with set operators "|" for union, "&" for intersection, and "-" for set difference.
        For example "trA-trB" would generate a Series with all domains of transcripts in "trA", but not in "trB".
    :param naming: Define how domains are named, either by "id" or by "name".
    :parm overlap_only: If set "True", only domains overlapping the region from column "start" to "end" are considered.
        If set "False", all domains of the transcripts are considered.
    :param insert_after: Define column after which the domains are inserted into the table, either by column name or index.
        By default, domain columns returned as seperate DataFrame.
    :param **filter_kwargs: additional keywords are passed to Gene.filter_transcripts, to restrict the transcripts to be considered.'''

    # set operators:
    # set union: |
    # set difference: -
    # set intersection: &

    # check arguments
    assert naming in ('id', 'name'), 'naming must be either "id" or "name".'
    label_idx = 0 if naming == 'id' else 1
    assert id_col in table, f'Missing id column "{id_col}" in table.'
    # check the "modes": can they be evaluated?
    tr_cols = {tr_col for mode in modes for tr_col in compile(mode, "<string>", "eval").co_names}
    for mode in modes:  # only set operations allowed, and should return a set
        assert isinstance(eval(mode, {tr_col: set() for tr_col in tr_cols}), set), f'{mode} does not return a set'
    # Do they contain only table names?
    missing = [c for c in tr_cols if c not in table.columns]
    assert len(missing) == 0, f'Missing transcript id columns in table: {", ".join(missing)}.'
    if insert_after is not None:
        if isinstance(insert_after, str):
            assert insert_after in table.columns, 'cannot find column "{insert_after}" in table.'
            insert_after = table.columns.get_loc(insert_after)
        assert isinstance(insert_after, int), 'insert_after must be a column name or column index'

    domain_rows = {}
    for idx, row in table.iterrows():
        domains = []
        g = transcriptome[row[id_col]]
        if filter_kwargs:
            valid_transcripts = set(g.filter_transcripts(**filter_kwargs))
        domain_sets = {}
        for tr_col in tr_cols:
            domain_sets[tr_col] = set()
            trids = set(row[tr_col]) & valid_transcripts if filter_kwargs else row[tr_col]
            for trid in trids:
                for dom in g.transcripts[trid].get('domain', {}).get(source, []):
                    if categories is not None and dom[2] not in categories:
                        continue
                    if overlap_only and not has_overlap(dom[4], (row.start, row.end)):
                        continue
                    domain_sets[tr_col].add(str(dom[label_idx]))
        for mode in modes:
            # evaluate string in mode
            domains.append(eval(mode, domain_sets))
        domain_rows[idx] = domains
    domain_rows = pd.DataFrame.from_dict(domain_rows, orient='index', columns=[f'{mode} domains' for mode in modes])
    if insert_after is None:
        return domain_rows
    return pd.concat([table.iloc[:, :insert_after+1], domain_rows, table.iloc[:, insert_after+1:]], axis=1)


def import_hmmer_models(path, model_file="Pfam-A.hmm.gz", metadata_file="Pfam-A.hmm.dat.gz"):
    '''Import the hmmer model and metadata.

    This function imports the hmmer Pfam models from "Pfam-A.hmm.gz" and metadata from "Pfam-A.hmm.dat.gz",
    which are available for download on the interpro website, at "https://www.ebi.ac.uk/interpro/download/Pfam/".
    The models are needed for Transcriptome.add_hmmer_domains() function.

    :param path: The path where model and metadata files are located.
    :param model_file: The filename of the model file.
    :param model_file: The filename of the metadata file.'''

    metadata = parse_hmmer_metadata(f'{path}/{metadata_file}')
    with pyhmmer.plan7.HMMFile((f'{path}/{model_file}')) as hmm_file:
        hmm_list = list(hmm_file)
    return metadata, hmm_list


def get_hmmer_sequences(transcriptome, genome_fn, aa_alphabet, query=True, ref_query=False, region=None, min_coverage=None,
                        max_coverage=None, gois=None, progress_bar=False):
    '''Get protein sequences in binary hmmer format.'''
    tr_ids = {}
    if query:
        for g, trids, _ in transcriptome.iter_transcripts(genewise=True, query=query, region=region, min_coverage=min_coverage,
                                                          max_coverage=max_coverage, gois=gois, progress_bar=progress_bar):
            tr_ids.setdefault(g.id, [[], []])[0] = trids
    if ref_query:
        for g, trids, _ in transcriptome.iter_ref_transcripts(genewise=True, query=ref_query, region=region, gois=gois, progress_bar=progress_bar):
            tr_ids.setdefault(g.id, [[], []])[1] = trids

    sequences = []
    seq_ids = []
    with FastaFile(genome_fn) as genome_fh:
        for gid in tr_ids:
            g = transcriptome[gid]
            seqs = {}
            for source in range(2):
                for trid, seq in g.get_sequence(genome_fh, tr_ids[gid][source], protein=True, reference=source).items():
                    seqs.setdefault(seq, []).append((gid, source, trid))
            for seq, seqnames in seqs.items():
                # Hack: use "name" attribute to store an integer.
                # Must be string encoded since it is interpreted as 0 terminated string and thus truncated
                text_seq = pyhmmer.easel.TextSequence(sequence=seq, name=bytes(str(len(sequences)), 'utf-8'))
                sequences.append(text_seq.digitize(aa_alphabet))
                seq_ids.append(seqnames)
    return sequences, seq_ids

#  function of isoseq.Transcriptome


def add_hmmer_domains(self, domain_models, genome, query=True, ref_query=False, region=None,
                      min_coverage=None, max_coverage=None, gois=None, progress_bar=False):
    '''Align domains to protein sequences  using pyhmmer and add them to the transcript isoforms.

    :param domain_models: The domain models and metadata, imported by "isotools.domains.import_hmmer_models" function
    :param genome: Filename of genome fasta file, or Fasta
    :param query: Query string to select the isotools transcripts, or True/False to include/exclude all transcripts.
    :param ref_query: Query string to select the reference transcripts, or True/False to include/exclude all transcripts.
    :param region: The region to be considered. Either a string "chr:start-end", or a tuple (chr,start,end). Start and end is optional.
    :param min_coverage: The minimum coverage threshold. Transcripts with less reads in total are ignored.
    :param max_coverage: The maximum coverage threshold. Transcripts with more reads in total are ignored.
    :param progress_bar: Print progress bars.
    '''

    metadata, models = domain_models
    # 1) get the protein sequences
    pipeline = pyhmmer.plan7.Pipeline(models[0].alphabet)
    logging.info('extracting protein sequences...')
    sequences, seq_ids = get_hmmer_sequences(self, genome, models[0].alphabet, query, ref_query, region=region,
                                             min_coverage=min_coverage, max_coverage=max_coverage, gois=gois,  progress_bar=False)
    logging.info(f'found {len(sequences)} different protein sequences from {sum(len(idL) for idL in seq_ids)} coding transcripts.')

    # 2) align domain models to sequences
    logging.info(f'aligning {len(models)} hmmer domain models to protein sequences...')
    hits = {hmm.accession.decode(): pipeline.search_hmm(hmm, sequences) for hmm in tqdm(models, disable=not progress_bar, unit='domains')}

    # 3) sort domains by gene/source/transcript
    domains = {}
    for pfam_acc, hL in hits.items():
        infos = metadata.loc[pfam_acc]
        for h in hL:
            seq_nr = int(h.name.decode())
            for domL in h.domains:
                ali = domL.alignment
                transcript_pos = (ali.target_from*3, ali.target_to*3)
                domains.setdefault(seq_nr, []).append(
                    (pfam_acc, infos['ID'], infos['TP'], transcript_pos, h.score, h.pvalue))
                # print(f'{h.name}\t{hmm.name}:\t{ali.target_from}-{ali.target_to}') #which sequence?

    # 4) add domains to transcripts
    dom_count = [0, 0]
    tr_count = [0, 0]
    for seq_nr, domL in domains.items():
        for gid, reference, trid in seq_ids[seq_nr]:
            g = self[gid]
            tr = g.ref_transcripts[trid] if reference else g.transcripts[trid]
            # get the genomic position of the domain boundaries
            orf = sorted(g.find_transcript_positions(trid, tr.get('CDS', tr.get('ORF'))[:2], reference=reference))
            pos_map = genomic_position([p+orf[0] for dom in domL for p in dom[3]], tr['exons'], g.strand == '-')
            trdom = tuple((*dom[:4], (pos_map[dom[3][0]+orf[0]], pos_map[dom[3][1]+orf[0]]), *dom[4:]) for dom in domL)
            tr.setdefault('domain', {})['hmmer'] = trdom
            tr_count[reference] += 1
            dom_count[reference] += len(domL)
    logger.info(f'found domains at {dom_count[1]} loci for {tr_count[1]} reference transcripts ' +
                f'and at {dom_count[0]} loci for {tr_count[0]} long read transcripts.')


def add_annotation_domains(self, annotation, category, id_col='uniProtId', name_col='name', inframe=True, progress_bar=False):
    '''Annotate isoforms with protein domains from uniprot ucsc table files.

    This function adds protein domains and other protein annotation to the transcripts.
    Annotation tables can be retriefed from https://genome.ucsc.edu/cgi-bin/hgTables. Select
    group: = "Genes and Gene Predictions", track: "UniProt", and chose from the available tables (e.g. "domains").

    :param annotation: The file name of the table downloaded from the table browser, or a pandas dataframe with the content.
    :param category: A term describing the type of annotation in the table (e.g. "domains").
    :param id_col: The column of the table with the domain id.
    :param name_col: The column of the table with the label.
    :param inframe: If set True (default), only annotations starting in frame are added to the transcript.
    :param append: If set True, the annotation is added to existing annotation. This may lead to duplicate entries.
        By default, annotation of the same category is removed before annotation is added.
    :param progress_bar: If set True, the progress is depicted with a progress bar.'''

    domain_count = 0
    # clear domains of that category
    if isinstance(annotation, str):
        anno = pd.read_csv(annotation, sep='\t')
    elif isinstance(annotation, pd.DataFrame):
        anno = annotation
    else:
        raise ValueError('"annotation" should be file name or pandas.DataFrame object')
    anno = anno.rename({'#chrom': 'chrom'}, axis=1)
    not_found = [c for c in ['chrom', 'chromStart', 'chromEnd', 'chromStarts', 'blockSizes', name_col] if c not in anno.columns]
    assert len(not_found) == 0, f'did not find the following columns in the annotation table: {", ".join(not_found)}'
    for _, row in tqdm(anno.iterrows(), total=len(anno), disable=not progress_bar, unit='domains'):
        if row['chrom'] not in self.chromosomes:
            continue
        for g in self.iter_genes(region=(row['chrom'], row.chromStart, row.chromEnd)):
            block_starts, block_sizes = list(map(int, row.chromStarts.split(','))), list(map(int, row.blockSizes.split(',')))
            blocks = []
            for s, l in zip(block_starts, block_sizes):
                if not blocks or row.chromStart+s > blocks[-1][1]:
                    blocks.append([row.chromStart+s, row.chromStart+s+l])
                else:
                    blocks[-1][1] = row.chromStart+s+l

            for ref in range(2):
                transcripts = g.ref_transcripts if ref else g.transcripts
                if not transcripts:
                    continue
                sg = g.ref_segment_graph if ref else g.segment_graph
                trids = [trid for trid in sg.search_transcript(blocks, complete=False, include_ends=True)
                         if 'ORF' in transcripts[trid] or 'CDS' in transcripts[trid]]
                for trid in trids:
                    tr = transcripts[trid]
                    try:
                        orf_pos = sorted(g.find_transcript_positions(trid, tr.get('CDS', tr.get('ORF'))[:2], reference=ref))
                        domain_pos = sorted(g.find_transcript_positions(trid, (row.chromStart, row.chromEnd), reference=ref))
                    except TypeError:  # > not supported for None, None
                        continue
                    if not (orf_pos[0] <= domain_pos[0] and domain_pos[1] <= orf_pos[1]):  # check within ORF
                        continue
                    if inframe and (domain_pos[0]-orf_pos[0]) % 3 != 0:  # check inframe
                        continue
                    domain_count += 1
                    dom_pos = (domain_pos[0]-orf_pos[0], domain_pos[1]-orf_pos[0])
                    dom_vals = (row[id_col], row[name_col], category,  dom_pos, (row.chromStart, row.chromEnd))
                    # check if present already
                    if dom_vals not in tr.setdefault('domain', {}).setdefault('annotation', []):
                        tr['domain']['annotation'].append(dom_vals)

    logger.info(f'found domains at {domain_count} transcript loci')


def get_interpro_domains(seqs, email, baseUrl='http://www.ebi.ac.uk/Tools/services/rest/iprscan5', progress_bar=True, max_jobs=25, poll_time=5):
    '''Request domains from ebi interpro REST API.

    Returns a list of the json responses as received, one for each requested sequenced.'''
    # examples at https://raw.githubusercontent.com/ebi-wp/webservice-clients/master/python/iprscan5.py

    requestUrl = baseUrl + u'/run/'
    current_jobs = {}
    if isinstance(seqs, str):
        seqs = [seqs]
    domains = [None]*len(seqs)
    i = 0
    with tqdm(unit='proteins', disable=not progress_bar, total=len(seqs)) as pbar:
        while seqs or current_jobs:
            while seqs and len(current_jobs) < max_jobs:
                params = {u'email': email, u'sequence': seqs.pop()}
                resp = requests.post(requestUrl, data=params)
                # todo: error handling - what if request fails?
                current_jobs[resp.content.decode()] = i
                pbar.set_description(f'waiting for {len(current_jobs)} jobs')
                i += 1
            time.sleep(poll_time)
            done = set()
            for job_id, idx in current_jobs.items():
                url = baseUrl + u'/status/' + job_id
                resp = requests.get(url)
                if not resp.ok:  # todo: error handling: e.g. timeout error?
                    continue
                status = resp.content.decode()
                if status in ('PENDING', 'RUNNING'):
                    continue
                if status == 'FINISHED':
                    url = baseUrl + u'/result/' + job_id + '/json'
                    resp = requests.get(url)
                    if resp.ok:
                        domains[idx] = resp.json()['results']  # else?
                        pbar.update(1)
                    else:
                        domains[idx] = [{'status': 'FAILED', 'reason': 'resp_not_ok', 'jobid': job_id}]
                    done.add(job_id)
                elif status == 'FAILED':  # try again?
                    logger.warning(f'Failed to get response for sequence {idx}')
                    done.add(job_id)
                    domains[idx] = [{'status': 'FAILED', 'reason': 'job failed', 'jobid': job_id}]
            for job_id in done:
                current_jobs.pop(job_id)
            pbar.set_description(f'waiting for {len(current_jobs)} jobs')
    return domains

# method of isotools.Gene


def add_interpro_domains(self, genome, email, baseUrl='http://www.ebi.ac.uk/Tools/services/rest/iprscan5',  max_jobs=25, poll_time=5,
                         query=True, ref_query=False, min_coverage=None, max_coverage=None,  progress_bar=True):
    '''Add domains to gene by webrequests to ebi interpro REST API.

    This function adds protein domains from interpro to the transcripts.
    Note that these rquest may take around 60 seconds per sequence.
    Seqeunces of transcripts sharing the same coding sequence are requested only once.

    :param genome: Filename of genome fasta
    :param email: Valid email address, as required by the ebi web request. Not sure how it is used.
    :param base_url: URL of the api.
    :param max_jobs: Specify the maximum number of parallel requests.
    :param query: Query string to select the isotools transcripts, or True/False to include/exclude all transcripts.
    :param ref_query: Query string to select the reference transcripts, or True/False to include/exclude all transcripts.
    :param min_coverage: The minimum coverage threshold. Transcripts with less reads in total are ignored.
    :param max_coverage: The maximum coverage threshold. Transcripts with more reads in total are ignored.
    :param progress_bar: If set True, the progress is depicted with a progress bar.'''

    seqs = {}  # seq -> trids dict, to avoid requesting the same sequence
    if query:
        trids = self.filter_transcripts(query=query, min_coverage=min_coverage, max_coverage=max_coverage)
        for trid, seq in self.get_sequence(genome, trids, protein=True).items():
            seqs.setdefault(seq, {}).setdefault('isotools', []).append(trid)
    if ref_query:
        ref_trids = self.filter_ref_transcripts(query=ref_query)
        for trid, seq in self.get_sequence(genome, ref_trids, protein=True, reference=True).items():
            seqs.setdefault(seq, {}).setdefault('reference', []).append(trid)

    dom_results = get_interpro_domains(list(seqs.keys()), email, baseUrl, progress_bar, max_jobs, poll_time)
    for i, (dom,) in enumerate(dom_results):
        if 'matches' not in dom:
            logger.warning(f'no response for sequence of {list(seqs.values())[i]}')
            continue
        domL = []
        for m in dom['matches']:
            for loc in m['locations']:
                entry = m['signature'].get('entry')
                domL.append((str(m['signature']['accession']),  # short name
                             str(m['signature']['name']),
                             entry.get('type', "unknown") if entry else "unknown",  # type
                             (loc['start']*3, loc['end']*3),  # position
                             loc.get('hmmBounds')))  # completeness
                # todo: potentially add more relevant information here
        for reference in range(2):
            for trid in seqs[dom['sequence']].get('reference' if reference else 'isotools', []):
                tr = self.ref_transcripts[trid] if reference else self.transcripts[trid]
                orf = sorted(self.find_transcript_positions(trid, tr.get('CDS', tr.get('ORF'))[:2], reference=reference))
                pos_map = genomic_position([p+orf[0] for d in domL for p in d[2]], tr['exons'], self.strand == '-')
                trdom = tuple((*d[:3], (pos_map[d[2][0]+orf[0]], pos_map[d[2][1]+orf[0]]), *d[3:]) for d in domL)
                tr.setdefault('domain', {})['interpro'] = trdom
