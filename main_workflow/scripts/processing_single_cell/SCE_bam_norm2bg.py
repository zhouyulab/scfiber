import argparse
import sys
from pybicl.io import iterline, BedFile
from pybicl.operation import IvalTree, Cluster
import pysam
from collections import defaultdict


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bam", required=True)
    base_group.add_argument("--ref", type=str, dest="ref", metavar="ref.gene.bed6", required=True)
    base_group.add_argument("--size-factor", type=str, dest="size_factor", metavar="size_factor", required=True)
    base_group.add_argument("--sample", type=str, dest="sample", metavar="sample", required=True)
    base_group.add_argument("--rep", type=str, dest="rep", metavar="rep", required=True)
    base_group.add_argument("--max-intron-len", type=int, dest="max_intron_len", metavar="max_intron_len",
                            required=False, default=1000000)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="prefix", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_ref = args.ref
    f_size_factor = args.size_factor
    select_sample = args.sample
    select_rep = args.rep
    max_intron_len = args.max_intron_len

    print("Loading reference gene info ...")
    ref_bed = BedFile(f_ref, "r")
    ref_ival = IvalTree(False)
    ref_ival.load_bed_file(ref_bed)

    print("Loading size factor info ...")
    size_factor_dict = dict()
    for indx, line in enumerate(iterline(f_size_factor)):
        if indx == 0:
            continue
        sample, rep, barcode, size_factor = line.rstrip("\n").split("\t")
        if sample != select_sample:
            continue
        if rep != select_rep:
            continue
        barcode = barcode.split("-")[0]
        size_factor = float(size_factor)
        size_factor_dict[barcode] = size_factor

    print("Loading alignments ...")
    bam = pysam.AlignmentFile(f_in, "rb")
    blocks_ival_dict = {x: IvalTree(True) for x in size_factor_dict.keys()}
    blocks_clu = Cluster(True)
    for read in bam.fetch():
        if read.is_secondary:
            continue
        if read.reference_end is None:
            continue
        if read.reference_start is None:
            continue
        if (read.reference_end - read.reference_start - read.query_length) > max_intron_len:
            continue
        barcode = read.get_tag("CR")
        chrom = read.reference_name
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"
        for s, e in read.blocks:
            blocks_ival_dict[barcode].add((s, e), chrom, s, e, strand)
            blocks_clu.add(chrom, s, e, strand)

    print("Normalization ...")
    f = {
        "+": open(f_out + ".Forward.bedGraph", "w"),
        "-": open(f_out + ".Reverse.bedGraph", "w")
    }

    for chrom, clu_start, clu_end, strand in blocks_clu.iter_all():
        overlap_blocks = list()
        for barcode, blocks_ival in blocks_ival_dict.items():
            tmp_overlap_blocks = blocks_ival.find(chrom, clu_start, clu_end, strand)
            if not tmp_overlap_blocks:
                continue
            size_factor = size_factor_dict[barcode]
            overlap_blocks += [(s, e, size_factor) for (s, e) in tmp_overlap_blocks]
        overlap_blocks = sorted(overlap_blocks, key=lambda x: [x[0], x[1]])
        boundaries = sorted({x[0] for x in overlap_blocks} | {x[1] for x in overlap_blocks})
        element_dict = {x: {"start": boundaries[x], "end": boundaries[x+1], "score": 0} for x in range(len(boundaries)-1)}
        boundary_indx = {x: boundaries.index(x) for x in boundaries}
        for s, e, size_factor in overlap_blocks:
            for indx in range(boundary_indx[s], boundary_indx[e]):
                element_dict[indx]["score"] += size_factor
        for val in element_dict.values():
            f[strand].write("{0}\t{1}\t{2}\t{3}\n".format(chrom, val["start"], val["end"], val["score"]))
    f["+"].close()
    f["-"].close()


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()