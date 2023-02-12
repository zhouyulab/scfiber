import argparse, sys, os
import pysam
from pybicl.io import iterline
from mpi4py import MPI
import math
from collections import defaultdict

class Peak(object):
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.barcode_set = defaultdict(set)

    def fetch_reads(self, bam):
        for read in bam.fetch(self.chrom, self.start, self.end):
            if read.mate_is_unmapped:
                continue
            if read.is_qcfail:
                continue
            if read.is_secondary:
                continue
            pcr_barcode = read.get_tag("BC")
            cell_barcode = read.get_tag("CR")
            self.barcode_set[cell_barcode].add(pcr_barcode)

    def count(self):
        count_dict = defaultdict(int)
        for cell_barcode, pcr_barcode_set in self.barcode_set.items():
            count_dict[cell_barcode] = len(pcr_barcode_set)
        return count_dict


def split_list(listTemp, n):
    length = len(listTemp)
    for i in range(n):
        yield listTemp[math.floor(i / n * length):math.floor((i + 1) / n * length)]

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bam", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="count.tsv", required=True)
    base_group.add_argument("--peak", type=str, dest="peak", metavar="peak.bed", required=True)
    return parser.parse_args(args)


def main(args):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rank_size = comm.Get_size()

    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_peak = args.peak

    if rank == 0:
        peak_li = list()
        for line in iterline(f_peak):
            data = line.rstrip("\n").split("\t")
            peak_li.append((data[0], int(data[1]), int(data[2])))
        split_peak_li = list(split_list(peak_li, rank_size))
    else:
        split_peak_li = None
    local_peak_li = comm.scatter(split_peak_li, root=0)
    bam = pysam.AlignmentFile(f_in, "rb")

    with open(f_out + ".tmp.{0}.tsv".format(rank), "w") as f:
        for key in local_peak_li:
            chrom, start, end = key
            peak = Peak(chrom, start, end)
            peak.fetch_reads(bam)
            cnt_str = ",".join(["{0}:{1}".format(bcd, cnt) for (bcd, cnt) in peak.count().items()])
            peak_str = "{0}_{1}_{2}".format(chrom, start, end)
            f.write("{0}\t{1}\n".format(peak_str, cnt_str))

    all_umi_cnt = comm.gather("done", root=0)

    if rank == 0:
        cell_barcodes = set()
        peaks = set()
        all_umi_peak_dict = dict()
        for rank_indx in range(rank_size):
            f_tmp = f_out + ".tmp.{0}.tsv".format(rank_indx)
            for line in iterline(f_tmp):
                peak, cnt_str = line.rstrip("\n").split("\t")
                peaks.add(peak)
                cnt_dict = defaultdict(int)
                if cnt_str:
                    for barcode_str in cnt_str.split(","):
                        bcd, cnt = barcode_str.split(":")
                        cnt_dict[bcd] = int(cnt)
                cell_barcodes |= set(cnt_dict.keys())
                all_umi_peak_dict[peak] = cnt_dict
            os.remove(f_tmp)

        cell_barcodes = sorted(cell_barcodes)
        with open(f_out, "w") as f:
            header = ["Peak"] + cell_barcodes
            f.write("\t".join(header) + "\n")
            for peak in sorted(peaks):
                peak_cnt_dict = all_umi_peak_dict[peak]
                cnt_li = [str(peak_cnt_dict[b]) for b in cell_barcodes]
                data = [peak] + cnt_li
                f.write("\t".join(data) + "\n")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
