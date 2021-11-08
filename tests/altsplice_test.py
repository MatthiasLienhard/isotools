# import isotools
# import pytest
# from isotools._transcriptome_filter import SPLICE_CATEGORY


def test_subcategory(example_gene):
    sg = example_gene.ref_segment_graph
    for novel in example_gene.transcripts:
        alt_splice = sg.get_alternative_splicing(novel['exons'])
        assert novel['transcript_name'] in alt_splice[1]
