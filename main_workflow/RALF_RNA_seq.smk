import os
curDir = os.getcwd()
include: "MetaInfo.sm.rules"

RALF_RNA_SEQ_BASE = os.path.join("analysis_v3", "RALF_RNA_seq")
RALF_SAMPLEs = ["RALF", "WT"]
RALF_REPs = ["rep1", "rep2"]

rule rm_adaptor:
    input:
        R1 = os.path.join("data", "RALF_RNA_seq", "{sample}.{rep}.R1.fq.gz"),
        R2 = os.path.join("data", "RALF_RNA_seq", "{sample}.{rep}.R2.fq.gz"),
    output:
        R1 = os.path.join(RALF_RNA_SEQ_BASE, "trim", "{sample}.{rep}.R1.fastq.gz"),
        R2 = os.path.join(RALF_RNA_SEQ_BASE, "trim", "{sample}.{rep}.R2.fastq.gz"),
    shell:
        """
cutadapt --cores 28 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 40 -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """

rule RNA_seq_mapping:
    input:
        trimed_R1 = rules.rm_adaptor.output.R1,
        trimed_R2 = rules.rm_adaptor.output.R2,
        gtf = REF_GTF,
    output:
        bam = os.path.join(RALF_RNA_SEQ_BASE, "bam", "{sample}.{rep}.sort.bam"),
    params:
        star_rRNA_indx = STAR_INDEX_RRNA,
        star_indx = STAR_INDX,
        mapping_dir = os.path.join(RALF_RNA_SEQ_BASE, "bam", "{sample}.{rep}"),
    shell:
        r"""
export LANG=en_EN
if [[ -e {params.mapping_dir}_rRNA ]]; then
    rm -rf {params.mapping_dir}_rRNA
fi
mkdir -p {params.mapping_dir}_rRNA
STAR --runMode alignReads --runThreadN 28 --genomeDir {params.star_rRNA_indx} \
    --sjdbOverhang 140 --outSAMmultNmax 1 --alignEndsType Local \
    --readFilesIn {input.trimed_R1} {input.trimed_R2} \
    --outFileNamePrefix {params.mapping_dir}_rRNA/ \
    --outReadsUnmapped Fastx --readFilesCommand gunzip -c

### Sort fastq code was from "https://www.biostars.org/p/59707/"
mkfifo {params.mapping_dir}_rRNA/tmp
awk 'NR%4==1{{n=$1}}NR%4==2{{s=$1}}NR%4==0{{print n,s,$1}}' {params.mapping_dir}_rRNA/Unmapped.out.mate1 | sort -S 2G > {params.mapping_dir}_rRNA/tmp &
awk 'NR%4==1{{n=$1}}NR%4==2{{s=$1}}NR%4==0{{print n,s,$1}}' {params.mapping_dir}_rRNA/Unmapped.out.mate2 | sort -S 2G | join -a1 -a2 {params.mapping_dir}_rRNA/tmp - | awk 'NF==5{{print $1"\n"$2"\n+\n"$3 >"{params.mapping_dir}_rRNA/R1.fq";print $1"\n"$4"\n+\n"$5 >"{params.mapping_dir}_rRNA/R2.fq"}}NF==3{{print $1"\n"$2"\n+\n"$3>"{params.mapping_dir}_rRNA/orphan.fq"}}'
rm {params.mapping_dir}_rRNA/tmp

if [[ -e {params.mapping_dir} ]]; then
    rm -rf {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runMode alignReads --runThreadN 28 --genomeDir {params.star_indx} \
--sjdbGTFfile {input.gtf} --sjdbOverhang 140 --outFilterMultimapNmax 5 --outSAMmultNmax 1 --alignEndsType Local \
--outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --alignIntronMax 100000 --alignMatesGapMax 100000 \
--readFilesIn {params.mapping_dir}_rRNA/R1.fq {params.mapping_dir}_rRNA/R2.fq \
--outFileNamePrefix {params.mapping_dir}/ --outReadsUnmapped Fastx
samtools view -@ 28 -q 20 -bS {params.mapping_dir}/Aligned.out.sam > {params.mapping_dir}/Aligned.out.bam
samtools sort -@ 28 {params.mapping_dir}/Aligned.out.bam -o {output.bam}
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam {params.mapping_dir}/Aligned.out.bam &
rm -r {params.mapping_dir}_rRNA
        """

rule rna_seq_rmDup:
    input:
        bam = rules.RNA_seq_mapping.output.bam,
    output:
        bam = os.path.join(RALF_RNA_SEQ_BASE, "rmDup", "{sample}.{rep}.rmDup.bam"),
        metrics = os.path.join(RALF_RNA_SEQ_BASE, "rmDup", "{sample}.{rep}.metrics"),
    shell:
        """
source activate py27
picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics}
samtools index {output.bam}
        """
        
rule rna_seq_bam2bw:
    input:
        bam = rules.rna_seq_rmDup.output.bam,
        merge_bed = REF_GENE_BED,
        size = CHROM_SIZE,
    output:
        flag = touch(os.path.join(RALF_RNA_SEQ_BASE, "bw", "{sample}.{rep}.flag")),
    params:
        total_sum = 1000000000,
        work_dir = os.path.join(RALF_RNA_SEQ_BASE, "bw", "{sample}.{rep}"),
    shell:
        """

rule_text=$(infer_experiment.py -i {input.bam} -r {input.merge_bed} | grep -e "[+-]")

indx=0
for word in ${{rule_text[@]}}; do
    let "indx = indx + 1"
    if [[ ${{indx}} = 6 ]]; then
        rule1=${{word:1:-2}}
    elif [[ ${{indx}} = 7 ]]; then
        score1=${{word}}
    elif [[ ${{indx}} = 13 ]]; then
        rule2=${{word:1:-2}}
    elif [[ ${{indx}} = 14 ]]; then
        score2=${{word}}
    fi
done

if [[ `awk "BEGIN{{print(${{score1}}>${{score2}})?"1":"0"}}"` = 1 ]]; then
    rule=${{rule1}}
    score=${{score1}}
    unrule_score=${{score2}}
else
    rule=${{rule2}}
    score=${{score2}}
    unrule_score=${{score1}}
fi

echo "file: {input.bam}  rule: ${{rule}}  score: ${{score}}"
if [[ `awk "BEGIN{{print(${{score}}-${{unrule_score}}>0.3)?"1":"0"}}"` = 0 ]]; then
    mode=1
else
    mode=2
fi

if [[ ${{mode}} = 1 ]]; then
    bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}
    wigToBigWig -clip {params.work_dir}.wig {input.size} {params.work_dir}.bw
    rm {params.work_dir}.wig
else
    bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir} -d "${{rule}}"
    wigToBigWig -clip {params.work_dir}.Forward.wig {input.size} {params.work_dir}.Forward.bw
    wigToBigWig -clip {params.work_dir}.Reverse.wig {input.size} {params.work_dir}.Reverse.bw
    rm {params.work_dir}.Forward.wig {params.work_dir}.Reverse.wig
fi
        """

rule RNA_feature_counts:
    input:
        bam = rules.rna_seq_rmDup.output.bam,
        gtf = REF_GTF,
    output:
        res = os.path.join(RALF_RNA_SEQ_BASE, "feature_count", "{sample}.{rep}.featurecount"),
    threads: 4
    shell:
        """
source activate py35
featureCounts -M -T {threads} -a {input.gtf} -o {output.res}.fc.out {input.bam}
cat {output.res}.fc.out | sed -n '3,$p' | awk '{{FS=OFS="\t"}}{{print $1, $7}}' > {output.res}
rm {output.res}.fc.out
        """

rule merge_RNA_feature_count_fpkm:
    input:
        flag = expand(rules.RNA_feature_counts.output.res, sample=RALF_SAMPLEs, rep=RALF_REPs),
    output:
        merge = os.path.join(RALF_RNA_SEQ_BASE, "merge_expr", "feature_counts.merge.tsv"),
    params:
        merge_featurecount = "Rscript scripts/RNA_seq/merge_featurecount.R",
    run:
        file_li = list()
        name_li = list()
        for sample in RALF_SAMPLEs:
            for rep in RALF_REPs:
                file_li.append(rules.RNA_feature_counts.output.res.format(sample=sample, rep=rep))
                name_li.append("{0}.{1}".format(sample, rep))
        cmd = "{merge_featurecount} -i {files} --name {name} -o {output}".format(
            merge_featurecount=params.merge_featurecount,
            files=" ".join(file_li), name=" ".join(name_li), output=output.merge
        )
        shell(cmd)

rule RALF_DEGs:
    input:
        fc = rules.merge_RNA_feature_count_fpkm.output.merge,
    output:
        res = os.path.join(RALF_RNA_SEQ_BASE, "DEGs", "RALF.DEGs.tsv"),
    params:
        RALF_DEGs = "Rscript scripts/RNA_seq/RALF_DEGs.R",
        out_dir = os.path.join(RALF_RNA_SEQ_BASE, "DEGs"),
        padj = 0.05,
        log2fc = 1,
    shell:
        """
{params.RALF_DEGs} {input.fc} {params.padj} {params.log2fc} {params.out_dir}
        """

rule RALF_DEGs_topGO:
    input:
        DE = rules.RALF_DEGs.output.res,
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RALF_RNA_SEQ_BASE, "DEGs", "topGO", "topGO.flag"),)
    params:
        work_dir = os.path.join(RALF_RNA_SEQ_BASE, "DEGs", "topGO"),
        run_topgo = "Rscript scripts/stat/run_topgo.R",
    shell:
        """
if [[ -e {params.work_dir} ]]; then
    rm -r {params.work_dir}
fi
mkdir -p {params.work_dir}/Up
mkdir -p {params.work_dir}/Down

cat {input} | awk '{{ FS=OFS="\t" }}{{ if($8=="Up"){{ print $7 }} }}' > {params.work_dir}/Up.gene.txt
cat {input} | awk '{{ FS=OFS="\t" }}{{ if($8=="Down"){{ print $7 }} }}' > {params.work_dir}/Down.gene.txt
{params.run_topgo} --gene-list {params.work_dir}/Up.gene.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/Up
{params.run_topgo} --gene-list {params.work_dir}/Down.gene.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/Down
        """

rule RALF_rMATS_turbo:
    input:
        WT_bam = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample="WT", rep=rep) for rep in RALF_REPs],
        RALF_bam = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample="RALF", rep=rep) for rep in RALF_REPs],
        gtf = REF_GTF,
    output:
        flag = touch(os.path.join(RALF_RNA_SEQ_BASE, "rMATS_turbo", "rmats.flag")),
    params:
        rmats = "python /home/wangdehe/downloads/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py",
        out_folder = os.path.join(RALF_RNA_SEQ_BASE, "rMATS_turbo"),
    threads: 28
    run:
        WT_cfg = os.path.join(params.out_folder, "WT_cfg.txt")
        with open(WT_cfg, "w") as f:
            f.write(",".join(input.WT_bam))

        RALF_cfg = os.path.join(params.out_folder, "RALF_cfg.txt")
        with open(RALF_cfg, "w") as f:
            f.write(",".join(input.RALF_bam))

        cmd="""
source activate py27
{rmats} --b1 {WT_cfg} --b2 {RALF_cfg} --nthread {threads} --gtf {gtf} --od {out_folder} -t paired --readLength 100
        """.format(
    rmats=params.rmats,
    WT_cfg=WT_cfg,
    RALF_cfg=RALF_cfg,
    gtf=input.gtf,
    out_folder=params.out_folder,
    threads=threads
)
        shell(cmd)

rule RALF_merge_bam:
    input:
        bams = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample=wildcards.sample, rep=rep) for rep in RALF_REPs],
    output:
        merged_bam = os.path.join(RALF_RNA_SEQ_BASE, "merged_bam", "{sample}.merged.bam"),
    shell:
        """
samtools merge -@ 50 -c -p {output.merged_bam}.tmp.bam {input.bams}
samtools sort -@ 50 -o {output.merged_bam} {output.merged_bam}.tmp.bam
samtools index {output.merged_bam}
rm {output.merged_bam}.tmp.bam
        """