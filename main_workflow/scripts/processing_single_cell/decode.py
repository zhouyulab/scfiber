import argparse, sys, os
import pysam
from pybicl.io import iterline


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bam", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output", required=True)
    base_group.add_argument("-b", "--barcode", type=str, dest="barcode", metavar="barcode.tsv", required=True)
    base_group.add_argument("-s", "--sample", type=str, dest="sample", metavar="sample", required=True)
    base_group.add_argument("-r", "--rep", type=str, dest="rep", metavar="rep", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_barcode = args.barcode
    select_sample = args.sample
    select_rep = args.rep

    barcode_clu_dict = dict()
    clu_set = set()
    for indx, line in enumerate(iterline(f_barcode)):
        if indx == 0:
            continue
        sample, rep, cluster, barcode = line.rstrip("\n").split("\t")
        if sample != select_sample:
            continue
        if rep != select_rep:
            continue
        if "-" in barcode:
            barcode = barcode.split("-")[0]
        barcode_clu_dict[barcode] = cluster
        clu_set.add(cluster)

    bam = pysam.AlignmentFile(f_in, "rb")
    clu_bam_out = {x: pysam.AlignmentFile(os.path.join(f_out, "{0}.bam".format(x)), "wb", header=bam.header) for x in clu_set}
    for read in bam.fetch():
        barcode = read.get_tag("CR")
        if barcode in barcode_clu_dict.keys():
            cluster = barcode_clu_dict[barcode]
            clu_bam_out[cluster].write(read)
            
    bam.close()
    for f in clu_bam_out.values():
        f.close()


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
