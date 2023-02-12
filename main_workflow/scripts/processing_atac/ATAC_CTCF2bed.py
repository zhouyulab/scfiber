import argparse, sys, os
from pybicl.io import iterline


class Hit(object):
    STRAND_DICT = {"n": "-", "p": "+"}

    def __init__(self, motif):
        self.motif = motif
        self.hit = list()

    def parse(self, line, f):
        tag, seq, peak, start, motif_len, strand, score = [x for x in map(lambda x: x.rstrip(";"), line.rstrip().split()) if x]
        start = int(start)
        motif_len = int(motif_len)
        score = float(score)
        peak_data = peak.split("_")
        chrom = "_".join(peak_data[:-2])
        peak_start = int(peak_data[-2])
        motif_start = peak_start + start
        motif_end = motif_start + motif_len
        data = [chrom, motif_start, motif_end, "{0}:{1}".format(peak, self.motif), score, self.STRAND_DICT[strand]]
        f.write("\t".join(list(map(str, data)))+"\n")


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.txt", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.bed6", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_out = args.output

    with open(f_out, "w") as f:
        begin = False
        for line in iterline(f_in):
            if line.startswith("AC"):
                assert not begin
                begin = True
                motif = Hit(line.rstrip().split()[-1])
            elif line.startswith("BS"):
                assert begin
                motif.parse(line, f)
            elif line.startswith("//"):
                assert begin
                begin = False
                del motif


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
