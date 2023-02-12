import argparse
import sys
from pybicl.io import BedFile, iterline
from pybicl.operation import IvalTree

def cnt_peak_base(peak_li, ref_ival):
    peak_base_cnt_dict = {
        "Intergenic": 0,
        "Promoter_1_2k": 0,
        "Promoter_0_1k": 0,
        "UTR5_exon": 0,
        "CDS_exon": 0,
        "UTR3_exon": 0,
        "Intron": 0,
        "TES_0_1k": 0,
        "TES_1_2k": 0
    }
    peak_cnt_dict = {
        "Intergenic": 0,
        "Promoter_1_2k": 0,
        "Promoter_0_1k": 0,
        "UTR5_exon": 0,
        "CDS_exon": 0,
        "UTR3_exon": 0,
        "Intron": 0,
        "TES_0_1k": 0,
        "TES_1_2k": 0
    }
    for peak in peak_li:
        peak_len = peak.chromEnd - peak.chromStart
        base_li = ["Intergenic"] * peak_len
        overlap_ref_li = list(ref_ival.find(peak.chrom, peak.chromStart-2000, peak.chromEnd+2000))
        flag_set = set()
        # Promoter_1_2k & TES_1_2k
        for rec in overlap_ref_li:
            if rec.strand == "+":
                promoter_start = rec.chromStart - 2000 - peak.chromStart
                promoter_end = rec.chromStart - 1000 - peak.chromStart
                tes_start = peak.chromEnd + 1000 - peak.chromStart
                tes_end = peak.chromEnd + 2000 - peak.chromStart
            else:
                tes_start = rec.chromStart - 2000 - peak.chromStart
                tes_end = rec.chromStart - 1000 - peak.chromStart
                promoter_start = peak.chromEnd + 1000 - peak.chromStart
                promoter_end = peak.chromEnd + 2000 - peak.chromStart
            promoter_start = min(max(0, promoter_start), peak_len)
            promoter_end = min(max(0, promoter_end), peak_len)
            tes_start = min(max(0, tes_start), peak_len)
            tes_end = min(max(0, tes_end), peak_len)
            promoter_len = promoter_end - promoter_start
            tes_len = tes_end - tes_start
            if tes_len:
                base_li = base_li[:tes_start] + ["TES_1_2k"] * tes_len + base_li[tes_end:]
                flag_set.add("TES_1_2k")
            if promoter_len:
                base_li = base_li[:promoter_start] + ["Promoter_1_2k"] * promoter_len + base_li[promoter_end:]
                flag_set.add("Promoter_1_2k")

        # Promoter_0_1k & TES_0_1k
        for rec in overlap_ref_li:
            if rec.strand == "+":
                promoter_start = rec.chromStart - 1000 - peak.chromStart
                promoter_end = rec.chromStart - peak.chromStart
                tes_start = peak.chromEnd - peak.chromStart
                tes_end = peak.chromEnd + 1000 - peak.chromStart
            else:
                tes_start = rec.chromStart - 1000 - peak.chromStart
                tes_end = rec.chromStart - peak.chromStart
                promoter_start = peak.chromEnd - peak.chromStart
                promoter_end = peak.chromEnd + 1000 - peak.chromStart
            promoter_start = min(max(0, promoter_start), peak_len)
            promoter_end = min(max(0, promoter_end), peak_len)
            tes_start = min(max(0, tes_start), peak_len)
            tes_end = min(max(0, tes_end), peak_len)
            promoter_len = promoter_end - promoter_start
            tes_len = tes_end - tes_start
            if tes_len:
                base_li = base_li[:tes_start] + ["TES_0_1k"] * tes_len + base_li[tes_end:]
                flag_set.add("TES_0_1k")
            if promoter_len:
                base_li = base_li[:promoter_start] + ["Promoter_0_1k"] * promoter_len + base_li[promoter_end:]
                flag_set.add("Promoter_0_1k")

        ## intron
        for rec in overlap_ref_li:
            for js in rec.iter_junction():
                js_start = js.chromStart - peak.chromStart
                js_end = js.chromEnd - peak.chromStart
                js_start = min(max(0, js_start), peak_len)
                js_end = min(max(0, js_end), peak_len)
                js_len = js_end - js_start
                if js_len:
                    base_li = base_li[:js_start] + ["Intron"] * js_len + base_li[js_end:]
                    flag_set.add("Intron")

        # UTR3
        for rec in overlap_ref_li:
            tmp = rec.to_utr3()
            if tmp is None:
                continue
            for exon in tmp.blocks:
                exon_start = exon.chromStart - peak.chromStart
                exon_end = exon.chromEnd - peak.chromStart
                exon_start = min(max(0, exon_start), peak_len)
                exon_end = min(max(0, exon_end), peak_len)
                exon_len = exon_end - exon_start
                if exon_len:
                    base_li = base_li[:exon_start] + ["UTR3_exon"] * exon_len + base_li[exon_end:]
                    flag_set.add("UTR3_exon")

        # UTR5
        for rec in overlap_ref_li:
            tmp = rec.to_utr5()
            if tmp is None:
                continue
            for exon in tmp.blocks:
                exon_start = exon.chromStart - peak.chromStart
                exon_end = exon.chromEnd - peak.chromStart
                exon_start = min(max(0, exon_start), peak_len)
                exon_end = min(max(0, exon_end), peak_len)
                exon_len = exon_end - exon_start
                if exon_len:
                    base_li = base_li[:exon_start] + ["UTR5_exon"] * exon_len + base_li[exon_end:]
                    flag_set.add("UTR5_exon")

        # CDS
        for rec in overlap_ref_li:
            tmp = rec.to_cds()
            if tmp is None:
                continue
            for exon in tmp.blocks:
                exon_start = exon.chromStart - peak.chromStart
                exon_end = exon.chromEnd - peak.chromStart
                exon_start = min(max(0, exon_start), peak_len)
                exon_end = min(max(0, exon_end), peak_len)
                exon_len = exon_end - exon_start
                if exon_len:
                    base_li = base_li[:exon_start] + ["CDS_exon"] * exon_len + base_li[exon_end:]
                    flag_set.add("CDS_exon")

        for base in base_li:
            peak_base_cnt_dict[base] += 1
        if flag_set:
            for flag in flag_set:
                peak_cnt_dict[flag] += 1
        else:
            peak_cnt_dict["Intergenic"] += 1
    return peak_base_cnt_dict, peak_cnt_dict


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="peak.bed6", required=True)
    base_group.add_argument("-r", "--ref", type=str, dest="ref", metavar="ref.bed12", required=True)
    base_group.add_argument("-s", "--size", type=str, dest="size", metavar="chrom.size", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="stat.tsv", required=True)
    return parser.parse_args(args)


def main(args):
    args = parse_args(args)
    f_in = args.input
    f_ref = args.ref
    f_size = args.size
    f_out = args.output

    ref_bed = BedFile(f_ref, "r")
    ref_ival = IvalTree(False)
    bg_ival = IvalTree(False)
    for rec in ref_bed.load("isoform"):
        ref_ival.add_record(rec)
        bg_ival.add(0, rec.chrom, rec.chromStart-2000, rec.chromEnd+2000, None)
    bg_cluster = bg_ival.trans2cluster()

    peak_bed = BedFile(f_in, "r")
    
    all_chrom_len = 0
    for line in iterline(f_size):
        all_chrom_len += int(line.rstrip("\n").split("\t")[1])
    bg_base_num_dict, _ = cnt_peak_base(bg_cluster.iter_all_rec(), ref_ival)
    bg_base_num_dict["Intergenic"] = 0
    bg_base_num_dict["Intergenic"] = all_chrom_len - sum(bg_base_num_dict.values())
    

    peak_base_num_dict, peak_cnt_dict = cnt_peak_base(peak_bed.load(), ref_ival)

    with open(f_out, "w") as f:
        header = ["Location", "PeakBaseNum", "AllBaseNum", "PeakNum"]
        f.write("\t".join(header)+"\n")
        for loc, num in sorted(peak_base_num_dict.items(), key=lambda x: x[1], reverse=True):
            data = [loc, str(num), str(bg_base_num_dict[loc]), str(peak_cnt_dict[loc])]
            f.write("\t".join(data) + "\n")


def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
