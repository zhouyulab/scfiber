import argparse, sys, os
from pybicl.io import iterline
from mpi4py import MPI
import math
from collections import defaultdict

def split_list(listTemp, n):
    length = len(listTemp)
    for i in range(n):
        yield listTemp[math.floor(i / n * length):math.floor((i + 1) / n * length)]

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.tsv", required=True)
    base_group.add_argument("--peak2gene", type=str, dest="peak2gene", metavar="peak2gene.tsv", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rank_size = comm.Get_size()

    args = parse_args(args)
    f_in = args.input
    f_peak2gene = args.peak2gene
    f_out = args.output

    if rank == 0:
        gene_peak_dict = defaultdict(list)
        for indx, line in enumerate(iterline(f_peak2gene)):
            if indx == 0:
                continue
            gid, peak = line.rstrip("\n").split("\t")
            gene_peak_dict[gid].append(peak)
    else:
        gene_peak_dict = None
    gene_peak_dict = comm.bcast(gene_peak_dict, root=0)

    tmp_file = f_out + ".{0}.tmp.tsv".format(rank)
    peak_name_li = None
    with open(tmp_file, "w") as f:
        for indx, line in enumerate(iterline(f_in)):
            if indx == 0:
                peak_name_li = line.rstrip("\n").split("\t")[3:]
            else:
                if indx % rank_size != rank:
                    continue
                data = line.rstrip("\n").split("\t")
                motif = data[0]
                rev_motif = data[1]
                cnt_li = list(map(int, data[3:]))
                gene_cnt_li = list()
                for gid, peak_li in sorted(gene_peak_dict.items()):
                    gene_cnt_li.append(sum([cnt_li[peak_name_li.index(x)] for x in peak_li]))
                f.write("\t".join(list(map(str, [motif, rev_motif]+gene_cnt_li)))+"\n")

    comm.gather("Done", root=0)
    if rank == 0:
        kmer_dict = dict()
        for tmp_rank in range(rank_size):
            tmp_file = f_out + ".{0}.tmp.tsv".format(tmp_rank)
            for line in iterline(tmp_file):
                motif = line.split("\t")[0]
                kmer_dict[motif] = line
            os.remove(tmp_file)

        with open(f_out, "w") as f:
            header = ["Motif", "RevMotif"] + sorted(gene_peak_dict.keys())
            f.write("\t".join(header) + "\n")
            for motif, line in sorted(kmer_dict.items()):
                f.write(line)


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
