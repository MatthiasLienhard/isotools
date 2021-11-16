# This script subsets the human genome fasta and annotation gff file as well as bam files to some regions of interest, as a test case
# run with:
#   python prepare_data.py --genome $hg38 --annotation $gff --out_prefix example
#   python prepare_data.py --genome $hg38 --annotation ${sa}.bam --out_prefix example_${sa}
import pysam
import Bio.bgzf
import argparse
import random


def main():
    def proportion(x):
        # https://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
        try:
            x = float(x)
        except ValueError:
            raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
        return x

    parser = argparse.ArgumentParser(prog='prepare_data', description='subset genomic data for testing')
    parser.add_argument('--annotation', metavar='<file.gff.gz>', help='specify reference annotation [requires tabix index]')
    parser.add_argument('--genome', metavar='<genome.fasta>', help='specify reference genome file [requires fai index]')
    parser.add_argument('--alignment', metavar='<alignment.bam>', help='specify alignment bam file [requires bai index]')
    parser.add_argument('--subsample_alignment', metavar='<fraction>', help='subsample the alignment file', type=proportion)
    parser.add_argument('--regions', metavar='<chr:start-end[,chr:start,end[,...]]>',
                        default='chr2:214700000-216200000,chr8:22000000-23000000', help='regions to subset')
    parser.add_argument('--out_prefix', metavar='</output/directory/prefix>', default='./example', help='specify output path and prefix')
    args = parser.parse_args()

    # parse regions
    regions = [(regstr[:regstr.find(':')], )+tuple(int(i) for i in regstr[regstr.find(':')+1:].split('-')) for regstr in args.regions.split(',')]
    # get subsets
    random.seed(42)  # make sure the same reads are picked if a fraction is specified
    if args.annotation:
        subset_annotation(args.annotation, regions, args.out_prefix+'.gff.gz')
    if args.genome:
        subset_genome(args.genome, regions, args.out_prefix+'.fa')
    if args.alignment:
        subset_alignment(args.alignment, regions, args.out_prefix+'.bam', proportion=args.subsample_alignment)

    return 0  # no error


def subset_genome(genome_fn, regions, out_fn="example_genome.fa"):
    fai = []
    offset = 0
    line_length = 50
    with open(out_fn, "w") as outfh:
        with pysam.FastaFile(genome_fn) as genome_fh:
            for reg in regions:
                outfh.write(f'>{reg[0]}_part\n')
                seq = (genome_fh.fetch(*reg))
                offset += len(reg[0])+7
                fai.append((f'{reg[0]}_part', len(seq), offset, line_length, line_length+1))
                for line in (seq[i:i+line_length] for i in range(0, len(seq), line_length)):
                    outfh.write(line+'\n')
                    offset += len(line)+1
    with open(out_fn+'.fai', "w") as outfh:
        for idx, reg in zip(fai, regions):
            outfh.write('\t'.join(str(v) for v in idx)+'\n')
            print(f'extracted {idx[1]} bases from {reg[0]}:{reg[1]}-{reg[2]}')


def subset_alignment(bam_fn, regions, out_fn="example.bam", proportion=None):
    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': end-start, 'SN': chrom+'_part'} for chrom, start, end in regions]}
    with pysam.AlignmentFile(out_fn, "wb", header=header) as out_fh:
        with pysam.AlignmentFile(bam_fn, "rb") as align:
            for new_ref_id, reg in enumerate(regions):
                n_reads = 0
                for read in align.fetch(*reg):
                    if proportion is not None and random.random() > proportion:
                        continue  # todo: this does not consider read pairs.
                    if read.reference_start < reg[1] or read.reference_end > reg[2]:
                        continue
                    if read.is_paired:
                        if read.next_reference_id != read.reference_id or read.next_reference_start < reg[1] or read.next_reference_end > reg[2]:
                            continue
                        read.next_reference_start -= reg[1]
                        read.next_reference_id = new_ref_id
                    read.reference_start -= reg[1]
                    read.reference_id = new_ref_id
                    out_fh.write(read)
                    n_reads += 1
                print(f'extracted {n_reads} reads from {reg[0]}:{reg[1]}-{reg[2]}')
    pysam.index(out_fn)
    # avoid pylint complaints: https://github.com/pysam-developers/pysam/issues/819


def subset_annotation(gff_fn, regions, out_fn="example_annotation.gff.gz"):
    gff = pysam.TabixFile(gff_fn)
    out_str = ''
    for reg in regions:
        genes = set()
        for line in gff.fetch(*reg):
            ls = line.split(sep="\t")
            if ls[2] == 'gene' and int(ls[3]) > reg[1] and int(ls[4]) < reg[2]:
                info = dict([pair.split('=', 1) for pair in ls[8].rstrip(';').split(";")])
                genes.add(info['gene_id'])
        print(f'extracted {len(genes)} genes from {reg[0]}:{reg[1]}-{reg[2]}')
        for line in gff.fetch(*reg):
            ls = line.split(sep="\t")
            ls[0] += '_part'
            ls[3] = str(int(ls[3])-reg[1])
            ls[4] = str(int(ls[4])-reg[1])
            info = dict([pair.split('=', 1) for pair in ls[8].rstrip(';').split(";")])
            if info.get('gene_id', None) in genes:
                out_str += ('\t'.join(ls) + '\n')
    with Bio.bgzf.BgzfWriter(out_fn, "wb") as outfh:
        outfh.write(out_str)
    _ = pysam.tabix_index(out_fn, preset='gff', force=True)


if __name__ == '__main__':
    exit(main())
