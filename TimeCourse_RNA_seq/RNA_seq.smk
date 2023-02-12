import os
curDir = os.getcwd()
include: "MetaInfo.smk"

SAMPLEs = ["WT", "FL"]
REPs = ["rep1", "rep2"]
TPs = ["n48h", "n36h", "n24h", "n12h", "0h", "12h", "24h", "36h", "48h"]

RNA_SEQ_BASE = os.path.join("analysis", "RNA_seq_TP")

rule RNA_cutadapt:
    input:
        R1 = os.path.join("data", "RNA_seq", "{sample}_{tp}.{rep}.R1.fq.gz"),
        R2 = os.path.join("data", "RNA_seq", "{sample}_{tp}.{rep}.R2.fq.gz"),
    output:
        R1 = os.path.join(RNA_SEQ_BASE, "trim", "{sample}", "{tp}.{rep}.R1.fastq.gz"),
        R2 = os.path.join(RNA_SEQ_BASE, "trim", "{sample}", "{tp}.{rep}.R2.fastq.gz"),
    threads: 28
    shell:
        """
cutadapt --cores 28 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 40 -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """

rule RNA_seq_mapping:
    input:
        trimed_R1 = rules.RNA_cutadapt.output.R1,
        trimed_R2 = rules.RNA_cutadapt.output.R2,
        gtf = REF_GTF,
    output:
        bam = os.path.join(RNA_SEQ_BASE, "bam", "{sample}", "{tp}.{rep}.sort.bam"),
    params:
        star_rRNA_indx = STAR_INDEX_RRNA,
        star_indx = STAR_INDX,
        mapping_dir = os.path.join(RNA_SEQ_BASE, "bam", "{sample}", "{tp}.{rep}"),
    threads: 28
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
        bam = os.path.join(RNA_SEQ_BASE, "rmDup", "{sample}", "{tp}.{rep}.rmDup.bam"),
        metrics = os.path.join(RNA_SEQ_BASE, "rmDup", "{sample}", "{tp}.{rep}.metrics"),
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
        chrom_size = CHROM_SIZE,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "bw", "{sample}", "{tp}.{rep}.flag")),
    params:
        total_sum = 1000000000,
        work_dir = os.path.join(RNA_SEQ_BASE, "bw", "{sample}", "{tp}.{rep}"),
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
    bam2wig.py -t {params.total_sum} -s {input.chrom_size} -i {input.bam} -o {params.work_dir}
    wigToBigWig -clip {params.work_dir}.wig {input.chrom_size} {params.work_dir}.bw
    rm {params.work_dir}.wig
else
    bam2wig.py -t {params.total_sum} -s {input.chrom_size} -i {input.bam} -o {params.work_dir} -d "${{rule}}"
    wigToBigWig -clip {params.work_dir}.Forward.wig {input.chrom_size} {params.work_dir}.Forward.bw
    wigToBigWig -clip {params.work_dir}.Reverse.wig {input.chrom_size} {params.work_dir}.Reverse.bw
    rm {params.work_dir}.Forward.wig {params.work_dir}.Reverse.wig
fi
        """

rule rna_seq_stringtie:
    input:
        bam = rules.rna_seq_rmDup.output.bam,
        merge_bed = REF_GENE_BED,
        gtf = REF_GTF,
    output:
        gtf = os.path.join(RNA_SEQ_BASE, "stringtie", "{sample}", "{tp}.{rep}.gtf"),
        tab = os.path.join(RNA_SEQ_BASE, "stringtie", "{sample}", "{tp}.{rep}.tab"),
    threads: 4
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
    rule_parm=""
else
    if [[ ${{rule}} = "++,--" || ${{rule}} == "1++,1--,2+-,2-+" ]]; then
        rule_parm="--rf"
    else
        rule_parm="--fr"
    fi
fi

stringtie -p 4 ${{rule_parm}} -e -G {input.gtf} -A {output.tab} -o {output.gtf} {input.bam}
        """

rule merge_stringtie_fpkm:
    input:
        flag = expand(rules.rna_seq_stringtie.output.gtf, sample=SAMPLEs, rep=REPs, tp=TPs),
    output:
        merge = os.path.join(RNA_SEQ_BASE, "merge_expr", "stringtie.merge.tsv"),
    params:
        merge_stringtie_fpkm = "Rscript scripts/RNA_seq/merge_stringtie_fpkm.R",
    run:
        file_li = list()
        name_li = list()
        for sample in SAMPLEs:
            for tp in TPs:
                for rep in REPs:
                    file_li.append(rules.rna_seq_stringtie.output.tab.format(sample=sample, tp=tp, rep=rep))
                    name_li.append("{0}.{1}.{2}".format(sample, tp, rep))
        cmd = "{merge_stringtie_fpkm} -i {files} --label {name} -o {output}".format(
            merge_stringtie_fpkm=params.merge_stringtie_fpkm,
            files=" ".join(file_li), name=" ".join(name_li), output=output.merge
        )
        shell(cmd)

rule RNA_feature_counts:
    input:
        bam = rules.rna_seq_rmDup.output.bam,
        gtf = REF_GTF,
    output:
        res = os.path.join(RNA_SEQ_BASE, "feature_count", "{sample}", "{tp}.{rep}.featurecount"),
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
        flag = expand(rules.RNA_feature_counts.output.res, sample=SAMPLEs, rep=REPs, tp=TPs),
    output:
        merge = os.path.join(RNA_SEQ_BASE, "merge_expr", "feature_counts.merge.tsv"),
    params:
        merge_featurecount = "Rscript scripts/RNA_seq/merge_featurecount.R",
    run:
        file_li = list()
        name_li = list()
        for sample in SAMPLEs:
            for tp in TPs:
                for rep in REPs:
                    file_li.append(rules.RNA_feature_counts.output.res.format(sample=sample, tp=tp, rep=rep))
                    name_li.append("{0}.{1}.{2}".format(sample, tp, rep))
        cmd = "{merge_featurecount} -i {files} --label {name} -o {output}".format(
            merge_featurecount=params.merge_featurecount,
            files=" ".join(file_li), name=" ".join(name_li), output=output.merge
        )
        shell(cmd)

rule find_circadian_gene:
    input:
        cnt = rules.merge_RNA_feature_count_fpkm.output.merge,
        interest_gene = os.path.join("data", "InterestGeneID.tsv"),
    output:
        circadian_group = os.path.join(RNA_SEQ_BASE, "circadian_gene", "CircadianGene.Group.tsv"),
        all_circadian_gene = os.path.join(RNA_SEQ_BASE, "circadian_gene", "AllCircadianGene.txt"),
        cpm = os.path.join(RNA_SEQ_BASE, "circadian_gene", "CPM.tsv"),
        cpm_merge = os.path.join(RNA_SEQ_BASE, "circadian_gene", "CPM.merge.tsv"),
    params:
        find_circadian_gene = "Rscript scripts/RNA_seq/find_circadian_gene.R",
        padj_cutoff = 0.01,
        lg2fc_cutoff = 0.5,
        out_dir = os.path.join(RNA_SEQ_BASE, "circadian_gene"),
    shell:
        """
{params.find_circadian_gene} {input.cnt} {input.interest_gene} {params.padj_cutoff} {params.lg2fc_cutoff} {params.out_dir}
        """

rule circadian_gene_term:
    input:
        circadian_group = rules.find_circadian_gene.output.circadian_group,
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "circadian_gene", "topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "circadian_gene", "topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
group=(
SignificantCircadian_Light
SignificantCircadian_Night
WeakCircadian_Light
WeakCircadian_Night
)
for group_indx in ${{group[@]}}; do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.circadian_group} | awk '$2=="'${{group_indx}}'"' | cut -f 1 > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule sig_circadian_gene_term:
    input:
        circadian_group = os.path.join(RNA_SEQ_BASE, "circadian_gene", "SigCircadianGene.Group.tsv"),
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "circadian_gene", "sig_topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "circadian_gene", "sig_topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
group=(
1
2
3
4
)
for group_indx in ${{group[@]}}; do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.circadian_group} | awk '$12=="'${{group_indx}}'"' | cut -f 1 > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule fetch_circadian_gene_promoter_seq:
    input:
        circadian_group = rules.find_circadian_gene.output.circadian_group,
        atac_peak = ATAC_PEAK,
        promoter_bed = PROMOTER_BED,
        genome_2bit = GENOME_2bit,
    output:
        circadian_group_gene_promoter_atac_bed = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "group{group_indx}.bed"),
        circadian_group_gene_promoter_atac_fa = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "group{group_indx}.fa"),
    params:
        circadian_group_gene_li = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "group{group_indx}.gene_li.txt"),
        circadian_group_gene_promoter_li = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "group{group_indx}.promoter.bed"),
    shell:
        """
cat {input.circadian_group} | awk '$2=={wildcards.group_indx}' | cut -f 1 > {params.circadian_group_gene_li}
cat {input.promoter_bed} | grep -f {params.circadian_group_gene_li} | sort -k 1,1 -k 2,2n > {params.circadian_group_gene_promoter_li}
bedtools intersect -wa -a {input.atac_peak} -b {params.circadian_group_gene_promoter_li} > {output.circadian_group_gene_promoter_atac_bed}
twoBitToFa -bed={output.circadian_group_gene_promoter_atac_bed} {input.genome_2bit} {output.circadian_group_gene_promoter_atac_fa}
        """

rule fetch_background_promoter_seq:
    input:
        atac_peak = ATAC_PEAK,
        promoter_bed = PROMOTER_BED,
        genome_2bit = GENOME_2bit,
    output:
        bg_promoter_atac_bed = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "background.bed"),
        bg_promoter_atac_fa = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "sequence", "background.fa"),
    shell:
        """
bedtools intersect -wa -a {input.atac_peak} -b {input.promoter_bed} > {output.bg_promoter_atac_bed}
twoBitToFa -bed={output.bg_promoter_atac_bed} {input.genome_2bit} {output.bg_promoter_atac_fa}
        """

rule circadian_gene_promoter_seq_homer:
    input:
        circadian_peak = rules.fetch_circadian_gene_promoter_seq.output.circadian_group_gene_promoter_atac_bed,
        bg_peak = rules.fetch_background_promoter_seq.output.bg_promoter_atac_bed,
        genome = GENOME_FA,
    output:
        motif = touch(os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "homer", "group{group_indx}.out", "done.flag")),
    params:
        out_dir = os.path.join(RNA_SEQ_BASE, "circadian_gene", "promoter_ATAC", "homer", "group{group_indx}.out")
    shell:
        """
findMotifsGenome.pl {input.circadian_peak} {input.genome} {params.out_dir} -bg {input.bg_peak}
        """

rule RNA_TP_rMATS_turbo:
    input:
        last_time_bam = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample=wildcards.sample, tp=wildcards.last_time, rep=rep) for rep in REPs],
        this_time_bam = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample=wildcards.sample, tp=wildcards.this_time, rep=rep) for rep in REPs],
        gtf = REF_GTF,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "rMATS_turbo", "{sample}", "{last_time}_vs_{this_time}", "rmats.flag")),
    params:
        rmats = "python /home/wangdehe/downloads/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py",
        out_folder = os.path.join(RNA_SEQ_BASE, "rMATS_turbo", "{sample}", "{last_time}_vs_{this_time}"),
    threads: 8
    run:
        last_time_cfg = os.path.join(params.out_folder, "last_time_cfg.txt")
        with open(last_time_cfg, "w") as f:
            f.write(",".join(input.last_time_bam))

        this_time_cfg = os.path.join(params.out_folder, "this_time_cfg.txt")
        with open(this_time_cfg, "w") as f:
            f.write(",".join(input.this_time_bam))

        cmd="""
source activate py27
{rmats} --b1 {last_time} --b2 {this_time} --nthread {threads} --gtf {gtf} --od {out_folder} -t paired --readLength 100
        """.format(
    rmats=params.rmats,
    last_time=last_time_cfg,
    this_time=this_time_cfg,
    gtf=input.gtf,
    out_folder=params.out_folder,
    threads=threads
)
        shell(cmd)

last_time_li = list()
this_time_li = list()
for indx in range(len(TPs)-1):
    last_time_li.append(TPs[indx])
    this_time_li.append(TPs[indx+1])

rule run_RNA_TP_rMATS_turbo:
    input:
        [rules.RNA_TP_rMATS_turbo.output.flag.format(sample=sample, last_time=TPs[last_indx], this_time=TPs[this_indx]) for sample in ["WT", "FL"] for this_indx in range(len(TPs)) for last_indx in range(this_indx)]

rule find_AS_circadian_gene:
    input:
        rMATS_flag = [rules.RNA_TP_rMATS_turbo.output.flag.format(sample=sample, last_time=last_time, this_time=this_time) for sample in SAMPLEs for (last_time, this_time) in zip(last_time_li, this_time_li)],
        interest_gene = os.path.join("data", "InterestGeneID.tsv"),
        cpm = rules.find_circadian_gene.output.cpm,
    output:
        AS_circadian_gene = os.path.join(RNA_SEQ_BASE, "AS_circadian_gene", "AllAsCircadianGene.txt"),
    params:
        find_AS_circadian_gene = "Rscript scripts/RNA_seq/find_AS_circadian_gene.R",
        FDR_cutoff = 0.05,
        delta_phi_cutoff = 0.2,
        read_num_cutoff = 20,
        rMATS_dir = os.path.join(RNA_SEQ_BASE, "rMATS_turbo"),
        out_dir = os.path.join(RNA_SEQ_BASE, "AS_circadian_gene"),
        tp_li = TPs,
        sample_li = SAMPLEs,
    shell:
        """
{params.find_AS_circadian_gene} -i {params.rMATS_dir} -t {params.tp_li} -s {params.sample_li} -o {params.out_dir} --interest-gene {input.interest_gene} --cpm {input.cpm} --FDR-cutoff {params.FDR_cutoff} --delta-phi-cutoff {params.delta_phi_cutoff} --read-num-cutoff {params.read_num_cutoff}
        """

rule AS_circadian_gene_term:
    input:
        AS_circadian_gene = rules.find_AS_circadian_gene.output.AS_circadian_gene,
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "AS_circadian_gene", "topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "AS_circadian_gene", "topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
if [[ -e {params.work_dir}/AllGene ]]; then
    rm -r {params.work_dir}/AllGene
fi
mkdir {params.work_dir}/AllGene
cat {input.AS_circadian_gene} | awk '{{print $2}}' | sort -u > {params.work_dir}/AllGene/gene_li.txt
{params.run_topgo} --gene-list {params.work_dir}/AllGene/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/AllGene

for group_indx in $(cat {input.AS_circadian_gene} | awk '{{print $1}}' | sort -u); do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.AS_circadian_gene} | awk '{{ if($1=="'${{group_indx}}'"){{print $2}} }}' | sort -u > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule find_time_course_gene:
    input:
        tp_fc = rules.merge_RNA_feature_count_fpkm.output.merge,
        tp_cpm = rules.find_circadian_gene.output.cpm_merge,
        tp_cpm_rep = rules.find_circadian_gene.output.cpm,
        interest_gene = os.path.join("data", "InterestGeneID.tsv"),
        tf_gene = os.path.join("data", "TF.gene_list.tsv"),
    output:
        course_gene = os.path.join(RNA_SEQ_BASE, "TimeCourse", "ImpulseDE2_enrich_gene.CPM.tsv")
    params:
        find_time_course_gene = "Rscript scripts/RNA_seq/find_time_course_gene.R",
        outdir = os.path.join(RNA_SEQ_BASE, "TimeCourse"),
    shell:
        """
{params.find_time_course_gene} --tp-fc {input.tp_fc} --tp-cpm {input.tp_cpm} --tp-cpm-rep {input.tp_cpm_rep} --interest-gene {input.interest_gene} --tf-gene {input.tf_gene} -o {params.outdir}
        """

rule time_course_gene_group_term:
    input:
        course_gene_gene = rules.find_time_course_gene.output.course_gene,
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "TimeCourse", "group_topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "TimeCourse", "group_topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
for group_indx in $(cat {input.course_gene_gene} | sed '1d' | awk '{{indx=NF-2; print $indx}}' | sort -u); do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.course_gene_gene} | awk '{{indx=NF-2; if($indx=="'${{group_indx}}'"){{print $1}} }}' | sort -u > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule time_course_gene_time_state_term:
    input:
        course_gene_gene = rules.find_time_course_gene.output.course_gene,
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "TimeCourse", "time_state_topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "TimeCourse", "time_state_topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
for group_indx in $(cat {input.course_gene_gene} | sed '1d' | awk '{{start_indx=NF-1; end_indx=NF; print $start_indx"_"$end_indx}}' | sort -u); do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.course_gene_gene} | awk '{{start_indx=NF-1; end_indx=NF; if($start_indx"_"$end_indx=="'${{group_indx}}'"){{print $1}} }}' | sort -u > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {input.topGO_bg} --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule TF_cor_net:
    input:
        tp_input = "analysis/RNA_seq_TP/merge_expr/feature_counts.merge.tsv",
        fc_tsv = [
    "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA2.merge.tsv",
    "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA4_part1.merge.tsv",
    "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA5.merge.tsv",
    "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA_web_part1.merge.tsv",
    "analysis/public_RNA_seq/merge_batch_expr/feature_counts.RNA_web_part2.merge.tsv",
    ],
        tf = os.path.join("data", "TF.gene_list.tsv"),
        interest_gene = os.path.join("data", "InterestGeneID.tsv"),
    output:
        res = expand(os.path.join(RNA_SEQ_BASE, "TfCorNet", "{stat}", "TF.CorNetGroup.tsv"), stat=["TP_pearson", "TP_spearman", "all_pearson", "all_spearman"])
    params:
        TF_cor_net = "Rscript scripts/RNA_seq/TF_cor_net.R",
        outdir = os.path.join(RNA_SEQ_BASE, "TfCorNet"),
    shell:
        """
{params.TF_cor_net} --tp-input {input.tp_input} --other-input {input.fc_tsv} --tf {input.tf} --interest-gene {input.interest_gene} -o {params.outdir}
        """

rule TF_cor_net_term:
    input:
        tp_cor_net = os.path.join(RNA_SEQ_BASE, "TfCorNet", "{stat}", "TF.CorNetGroup.tsv"),
        tf = os.path.join("data", "TF.gene_list.tsv"),
        topGO_map = TOP_GO_MAP,
        topGO_bg = TOP_GO_BG,
    output:
        flag = touch(os.path.join(RNA_SEQ_BASE, "TfCorNet", "{stat}", "topGO", "topGO.flag")),
    params:
        work_dir = os.path.join(RNA_SEQ_BASE, "TfCorNet", "{stat}", "topGO"),
        run_topgo = "Rscript scripts/utils/run_topgo.R",
    shell:
        """
cat {input.tf} | sed '1d' | awk '{{print $1}}' | sort -u > {params.work_dir}/background.txt
for group_indx in $(cat {input.tp_cor_net} | sed '1d' | awk '{{print $2}}' | sort -u); do
    if [[ -e {params.work_dir}/${{group_indx}} ]]; then
        rm -r {params.work_dir}/${{group_indx}}
    fi
    mkdir {params.work_dir}/${{group_indx}}
    cat {input.tp_cor_net} | awk '{{if($2=="'${{group_indx}}'"){{print $1}} }}' | sort -u > {params.work_dir}/${{group_indx}}/gene_li.txt
    {params.run_topgo} --gene-list {params.work_dir}/${{group_indx}}/gene_li.txt --background-list {params.work_dir}/background.txt --GO-db {input.topGO_map} --top-nodes 20 --output {params.work_dir}/${{group_indx}}
done
        """

rule run_TF_cor_net_term:
    input:
        expand(rules.TF_cor_net_term.output, stat=["TP_pearson", "TP_spearman", "all_pearson", "all_spearman"])

rule merge_bam:
    input:
        bams = lambda wildcards: [rules.rna_seq_rmDup.output.bam.format(sample=wildcards.sample, tp=wildcards.tp, rep=rep) for rep in REPs],
    output:
        merged_bam = os.path.join(RNA_SEQ_BASE, "merged_bam", "{sample}.{tp}.merged.bam"),
    shell:
        """
samtools merge -@ 50 -c -p {output.merged_bam}.tmp.bam {input.bams}
samtools sort -@ 50 -o {output.merged_bam} {output.merged_bam}.tmp.bam
samtools index {output.merged_bam}
rm {output.merged_bam}.tmp.bam
        """

rule merge_bam_tmp:
    input:
        expand(rules.merge_bam.output, sample=SAMPLEs, tp=TPs),

rule RNA_seq:
    input:
        expand(rules.rna_seq_bam2bw.output, sample=SAMPLEs, rep=REPs, tp=TPs),
        rules.merge_stringtie_fpkm.output,
        rules.merge_RNA_feature_count_fpkm.output,
        rules.find_circadian_gene.output,
        expand(rules.circadian_gene_promoter_seq_homer.output, group_indx=[1, 2, 3, 4]),
