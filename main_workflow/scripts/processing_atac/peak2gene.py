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

    f_peak = os.path.join(f_in, "peaks.bed")
    print("Loading peak data ...")
    peak_ival_tree = IvalTree(False)
    for indx, line in enumerate(iterline(f_peak)):
        chrom, start, end = line.rstrip("\n").split("\t")
        start = int(start)
        end = int(end)
        peak_ival_tree.add(indx, chrom, start, end, None)

    print("Loading barcode data ...")
    f_barcode = os.path.join(f_in, "barcodes.tsv")
    barcode_li = [x.rstrip("\n") for x in iterline(f_barcode)]

    print("Loading matrix data ...")
    f_mat = os.path.join(f_in, "matrix.mtx")
    cnt_mat = defaultdict(dict)
    for indx, line in enumerate(iterline(f_mat, remark="%")):
        if indx == 0:
            continue
        peak_indx, cell_indx, num = line.rstrip("\n").split(" ")
        peak_indx = int(peak_indx) - 1  # 1-base to 0-base
        cell_indx = int(cell_indx) - 1  # 1-base to 0-base
        num = int(num)
        cnt_mat[cell_indx][peak_indx] = num

    print("Peak to gene")
    ref_bed = BedFile(f_ref, "r")
    with open(f_out, "w") as f:
        header = "Gid\t" + "\t".join(barcode_li) + "\n"
        f.write(header)
        gene_num = 0
        for rec in ref_bed.load():
            gene_num += 1
            if gene_num % 100 == 0:
                print("writing {0} genes ...".format(gene_num))
            if rec.strand == "+":
                overlap_peak_indxs = set(peak_ival_tree.find(
                    rec.chrom, rec.chromStart - upstream_expand, rec.chromEnd + downstream_expand, None
                ))
            else:
                overlap_peak_indxs = set(peak_ival_tree.find(
                    rec.chrom, rec.chromStart - downstream_expand, rec.chromEnd + upstream_expand, None
                ))
            if overlap_peak_indxs:
                cnts = list()
                for cell_indx in range(len(barcode_li)):
                    peak_dict = cnt_mat[cell_indx]
                    if peak_dict:
                        select_peak_indx = overlap_peak_indxs & peak_dict.keys()
                        tmp_cnt = sum([peak_dict[x] for x in select_peak_indx])
                        cnts.append(tmp_cnt)
                    else:
                        cnts.append(0)
            else:
                cnts = [0 for _ in barcode_li]
            f.write("\t".join(list(map(str, [rec.name]+cnts)))+"\n")
    print("done")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
