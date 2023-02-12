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
    base_group.add_argument("--genome", type=str, dest="genome", metavar="genome.fa", required=True)
    base_group.add_argument("--size-factor", type=str, dest="size_factor", metavar="size_factor", required=True)
    base_group.add_argument("--sample", type=str, dest="sample", metavar="sample", required=True)
    base_group.add_argument("--rep", type=str, dest="rep", metavar="rep", required=True)
    base_group.add_argument("--max-intron-len", type=int, dest="max_intron_len", metavar="max_intron_len",
                            required=False, default=1000000)
    base_group.add_argument("--internal-polya-detected", type=int, dest="internal_polyA_detected", metavar="internal_polyA_detected",
                            required=False, default=8)
    base_group.add_argument("--internal-polya-cutoff", type=int, dest="internal_polyA_cutoff", metavar="internal_polyA_cutoff",
                            required=False, default=6)                       
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="prefix", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_ref = args.ref
    f_genome = genome = pysam.FastaFile(args.genome)
    f_size_factor = args.size_factor
    select_sample = args.sample
    select_rep = args.rep
    max_intron_len = args.max_intron_len
    internal_polyA_detected = args.internal_polyA_detected
    internal_polyA_cutoff = args.internal_polyA_cutoff

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
    polyA_dict = {"+": defaultdict(set), "-": defaultdict(set)}
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
            polyA = read.reference_start
            downstream_seq = genome.fetch(chrom, polyA, polyA+internal_polyA_detected).upper()
            if downstream_seq.count("A") >= internal_polyA_cutoff:
                continue
        else:
            strand = "+"
            polyA = read.reference_end - 1
            downstream_seq = genome.fetch(chrom, polyA-internal_polyA_detected, polyA).upper()
            if downstream_seq.count("T") >= internal_polyA_cutoff:
                continue
        polyA_dict[strand][(chrom, polyA)].add(barcode)

    print("Normalization ...")
    f = {
        "+": open(f_out + ".Forward.bedGraph", "w"),
        "-": open(f_out + ".Reverse.bedGraph", "w")
    }

    for strand, data in polyA_dict.items():
        for (chrom, polyA), barcodes in data.items():
            score = sum([size_factor_dict[barcode] for barcode in barcodes])
            f[strand].write("{0}\t{1}\t{2}\t{3}\n".format(chrom, polyA, polyA+1, score))
    f["+"].close()
    f["-"].close()


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
