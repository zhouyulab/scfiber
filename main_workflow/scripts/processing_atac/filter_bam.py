import argparse, sys, os
import pysam
from pybicl.io import iterline


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.bam", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bam", required=True)
    base_group.add_argument("-b", "--barcode", type=str, dest="barcode", metavar="barcode.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output
    f_barcode = args.barcode

    barcode_li = list()
    for indx, line in enumerate(iterline(f_barcode)):
        barcode = line.rstrip("\n")
        barcode_li.append(barcode)

    bam = pysam.AlignmentFile(f_in, "rb")
    bam_out = pysam.AlignmentFile(f_out, "wb", header=bam.header)
    for read in bam.fetch():
        barcode = read.get_tag("CR")
        if barcode in barcode_li:
            bam_out.write(read)
            
    bam.close()
    bam_out.close()


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
