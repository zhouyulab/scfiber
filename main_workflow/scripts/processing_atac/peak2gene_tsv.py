import argparse
import sys
import os
from pybicl.io import iterline, BedFile
from pybicl.operation import IvalTree
from collections import defaultdict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input_dir", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="ref.gene.bed6", required=True)
    base_group.add_argument("--upstream-expand", type=int, dest="upstream_expand", metavar="upstream_expand",
                            default=2000)
    base_group.add_argument("--downstream-expand", type=int, dest="downstream_expand", metavar="downstream_expand",
                            default=0)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_ref = args.ref
    upstream_expand = args.upstream_expand
    downstream_expand = args.downstream_expand

    print("Loading peak data ...")
    peak_ival_tree = IvalTree(False)
    for indx, line in enumerate(iterline(f_in)):
        chrom, start, end = line.rstrip("\n").split("\t")
        start = int(start)
        end = int(end)
        peak_ival_tree.add("{0}_{1}_{2}".format(chrom, start, end), chrom, start, end, None)

    print("Peak to gene")
    ref_bed = BedFile(f_ref, "r")
    with open(f_out, "w") as f:
        header = "Gid\tPeak\n"
        f.write(header)
        for rec in ref_bed.load():
            if rec.strand == "+":
                overlap_peaks = set(peak_ival_tree.find(
                    rec.chrom, rec.chromStart - upstream_expand, rec.chromEnd + downstream_expand, None
                ))
            else:
                overlap_peaks = set(peak_ival_tree.find(
                    rec.chrom, rec.chromStart - downstream_expand, rec.chromEnd + upstream_expand, None
                ))
            for peak in overlap_peaks:
                f.write("{0}\t{1}\n".format(rec.name, peak))

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
