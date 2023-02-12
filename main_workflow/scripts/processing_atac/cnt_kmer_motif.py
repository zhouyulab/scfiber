import argparse, sys, os
import pysam
from pybicl.io import iterline
from pybicl.utils import Interval
from mpi4py import MPI
import math
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

BASEs = ["A", "C", "G", "T"]

def enum_base(num, prefix=""):
    if num > 0:
        for b in BASEs:
            tmp_str = prefix + b
            for res in enum_base(num - 1, prefix=tmp_str):
                yield res
    else:
        yield prefix

def init_kmer_motif(kmer):
    motif_li = list()
    for kmer_motif in enum_base(kmer):
        rev_motif = Seq(kmer_motif, IUPAC.unambiguous_dna).reverse_complement()
        if str(rev_motif) in motif_li:
            continue
        motif_li.append(kmer_motif)
    return motif_li

def split_list(listTemp, n):
    length = len(listTemp)
    for i in range(n):
        yield listTemp[math.floor(i / n * length):math.floor((i + 1) / n * length)]

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--peak", type=str, dest="peak", metavar="peak.bed3", required=True)
    base_group.add_argument("-g", "--genome", type=str, dest="genome", metavar="genome.fa", required=True)
    base_group.add_argument("-k", "--kmer", type=int, dest="kmer", metavar="kmer", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="count.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rank_size = comm.Get_size()

    args = parse_args(args)
    f_peak = args.peak
    f_genome = args.genome
    f_out = args.output
    kmer = args.kmer

    if rank == 0:
        motif_li = init_kmer_motif(kmer)
        split_motif_li = list(split_list(motif_li, rank_size))
    else:
        split_motif_li = None
    local_motif_li = comm.scatter(split_motif_li, root=0)
    local_rev_motif_li = [str(Seq(x, IUPAC.unambiguous_dna).reverse_complement()) for x in local_motif_li]
    genome = pysam.FastaFile(f_genome)

    tmp_file = f_out + ".{0}.tmp.tsv".format(rank)
    peak_li = list()
    with open(tmp_file, "w") as f:
        f.write("\t".join(local_motif_li)+"\n")
        for line in iterline(f_peak):
            data = line.rstrip("\n").split("\t")
            chrom = data[0]
            start = int(data[1])
            end = int(data[2])
            peak_li.append("{0}_{1}_{2}".format(chrom, start, end))
            ival = Interval(chrom, start, end)
            peak_seq = ival.fetch_sequence(genome).upper()
            motif_cnt = list()
            for motif, rev_motif in zip(local_motif_li, local_rev_motif_li):
                tmp_cnt = peak_seq.count(motif)
                if motif != rev_motif:
                    tmp_cnt += peak_seq.count(rev_motif)
                motif_cnt.append(tmp_cnt)
            f.write("\t".join(list(map(str, motif_cnt)))+"\n")

    comm.gather("Done", root=0)
    if rank == 0:
        cnt_dict = defaultdict(list)
        for tmp_rank in range(rank_size):
            tmp_file = f_out + ".{0}.tmp.tsv".format(tmp_rank)
            tmp_motif_li = None
            for indx, line in enumerate(iterline(tmp_file)):
                data = line.rstrip("\n").split("\t")
                if indx == 0:
                    tmp_motif_li = data
                else:
                    for tmp_motif, tmp_cnt in zip(tmp_motif_li, data):
                        cnt_dict[tmp_motif].append(tmp_cnt)
            os.remove(tmp_file)

        with open(f_out, "w") as f:
            header = ["Motif", "RevMotif"] + peak_li
            f.write("\t".join(header) + "\n")
            peak_len = len(peak_li)
            for motif, cnts in sorted(cnt_dict.items()):
                assert len(cnts) == peak_len
                rev_motif = str(Seq(motif, IUPAC.unambiguous_dna).reverse_complement())
                data = [motif, rev_motif] + cnts
                f.write("\t".join(data)+"\n")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
