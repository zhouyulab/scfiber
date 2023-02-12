import argparse, sys
import pysam
from pybicl.operation import IvalTree
from pybicl.io import BedFile
from pybicl.utils import Interval
import re


COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "N": "N", "W": "W", "S": "S", "K": "M",
    "M": "K", "Y": "R", "R": "Y", "H": "D",
    "B": "V", "V": "B", "D": "H"
}


def rev_com(seq):
    return "".join([COMPLEMENT[x] for x in seq.upper()])[::-1]


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="refgene.bed6", required=True)
    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fa", required=True)
    base_group.add_argument("-p", "--peak", type=str, dest="peak", metavar="peak.bed6", required=True)
    base_group.add_argument("-m", "--motif", type=str, dest="motif", metavar="motif.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.gtf", required=True)
    base_group.add_argument("--upstream-expand", type=int, dest="up_expand", metavar="up_expand", default=1000)
    base_group.add_argument("--downstream-expand", type=int, dest="down_expand", metavar="down_expand", default=500)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_ref = args.ref
    f_peak = args.peak
    f_motif = args.motif
    f_genome = args.genome
    f_out = args.output
    up_expand = args.up_expand
    down_expand = args.down_expand

    ref_gene_bed = BedFile(f_ref, "r")
    genome = pysam.FastaFile(f_genome)
    peak_bed = BedFile(f_peak, "r")
    peak_ival = IvalTree(False)
    peak_ival.load_bed_file(peak_bed, "bed")
    motif_li = list()
    rev_motif_li = list()
    with open(f_motif, "r") as f:
        for line in f.readlines():
            motif = line.rstrip().upper()
            motif_li.append(motif)
            rev_motif_li.append(rev_com(motif))

    with open(f_out, "w") as f:
        header = [
            "Gid", "Chrom", "TSS", "Strand", "PeakStart", "PeakEnd",
            "Motif", "RevMotif", "CisMotif", "MotifStart", "MotifEnd"
        ]
        f.write("\t".join(header)+"\n")
        motif_cis_dict = {"+": "Motif", "-": "RevMotif"}
        rev_motif_cis_dict = {"-": "Motif", "+": "RevMotif"}
        for gene in ref_gene_bed.load("bed"):
            if gene.strand == "+":
                tss = gene.chromStart
                search_region = Interval(gene.chrom, tss - up_expand, tss + down_expand, None)
            else:
                tss = gene.chromEnd - 1
                search_region = Interval(gene.chrom, tss - down_expand, tss + up_expand, None)
            for peak in peak_ival.find_overlap(search_region):
                sequence = genome.fetch(peak.chrom, peak.chromStart, peak.chromEnd).upper()
                for motif, rev_motif in zip(motif_li, rev_motif_li):
                    for motif_hit in re.finditer(motif, sequence):
                        data = [
                            gene.name, gene.chrom, tss, gene.strand, peak.chromStart, peak.chromEnd,
                            motif, rev_motif, motif_cis_dict[gene.strand],
                            motif_hit.start() + peak.chromStart, motif_hit.end() + peak.chromStart
                        ]
                        f.write("\t".join(list(map(str, data)))+"\n")
                    for rev_hit in re.finditer(rev_motif, sequence):
                        data = [
                            gene.name, gene.chrom, tss, gene.strand, peak.chromStart, peak.chromEnd,
                            motif, rev_motif, rev_motif_cis_dict[gene.strand],
                            rev_hit.start() + peak.chromStart, rev_hit.end() + peak.chromStart
                        ]
                        f.write("\t".join(list(map(str, data))) + "\n")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()