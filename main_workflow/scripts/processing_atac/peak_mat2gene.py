import argparse
import sys
import os
from pybicl.io import iterline
import numpy as np
from pybicl.operation import IvalTree
from collections import defaultdict
from mpi4py import MPI
import math

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.tsv", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="ref.gene.bed6", required=True)
    base_group.add_argument("--upstream-expand", type=int, dest="upstream_expand", metavar="upstream_expand",
                            default=2000)
    base_group.add_argument("--downstream-expand", type=int, dest="downstream_expand", metavar="downstream_expand",
                            default=0)
    return parser.parse_args(args)

def split_list(listTemp, n):
    length = len(listTemp)
    for i in range(n):
        yield listTemp[math.floor(i / n * length):math.floor((i + 1) / n * length)]

def main(args):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rank_size = comm.Get_size()

    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_ref = args.ref
    upstream_expand = args.upstream_expand
    downstream_expand = args.downstream_expand

    if rank == 0:
        ref_gene_li = list(iterline(f_ref))
        split_ref = list(split_list(ref_gene_li, rank_size))
    else:
        split_ref = None
    local_ref_li = comm.scatter(split_ref, root=0)

    local_ref_chrom = set([x.split("\t")[0] for x in local_ref_li])

    print("Loading peak data ...")
    peak_ival_tree = IvalTree(False)
    peak_cnt_info = dict()
    barcods = None
    for indx, line in enumerate(iterline(f_in)):
        data = line.rstrip("\n").split("\t")
        key = data[0]
        value = data[1:]
        if indx == 0:
            barcods = value
        else:
            peak_name = key
            name_data = peak_name.split("_")
            chrom = "_".join(name_data[:-2])
            if chrom not in local_ref_chrom:
                continue
            start = int(name_data[-2])
            end = int(name_data[-1])
            cnts = np.array(list(map(int, value)))
            peak_ival_tree.add(peak_name, chrom, start, end, None)
            peak_cnt_info[peak_name] = cnts

    with open(f_out+".{0}.tmp.tsv".format(rank), "w") as f:
        for line in local_ref_li:
            rec_chrom, rec_start, rec_end, rec_name, _, rec_strand  = line.rstrip("\n").split("\t")
            
            if rec_strand == "+":
                overlap_peak_name = set(peak_ival_tree.find(
                    rec_chrom, int(rec_start) - upstream_expand, int(rec_end) + downstream_expand, None
                ))
            else:
                overlap_peak_name = set(peak_ival_tree.find(
                    rec_chrom, int(rec_start) - downstream_expand, int(rec_end) + upstream_expand, None
                ))
            cnts = np.array([0] * len(barcods))
            for peak_name in overlap_peak_name:
                cnts = cnts + peak_cnt_info[peak_name]
            f.write("\t".join(list(map(str, [rec_name]+list(cnts))))+"\n")
    
    flag = comm.gather("done", root=0)
    if rank == 0:
        with open(f_out, "w") as f:
            header = ["Gene"] + barcods
            f.write("\t".join(header)+"\n")
            for rank_indx in range(rank_size):
                f_tmp = f_out+".{0}.tmp.tsv".format(rank_indx)
                for line in iterline(f_tmp):
                    f.write(line)
                os.remove(f_tmp)


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
