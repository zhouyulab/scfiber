# scfiber

This repository is for the study of integrating scRNA-seq and scATAC-seq to systematically characterize the cells of the outer integument of ovules from wildtype and fuzzless/lintless cotton. We identified five cell populations including the fiber cell type, and uncovered fiber-specific circadian clock-controlled gene expression program in regulating fiber growth.

## Source Code Organization

The analysis workflows were implemented and managed using Snakemake.

- main_workflow: Main workflow
- SCAN_Fiber5DPA: Data processing of SCAN-seq and Fiber RNA-seq (5DPA)
- TimeCourse_RNA_seq: Data processing, DEGs and circadian gene analysis of time-course RNA-seq

## Reference

Wang et al. Cell-specific Clock Controlled Gene Expression Program Regulates 1 Rhythmic Fiber Cell Growth in Cotton. Accepted in Genome Biology. 2023.

## Developers

* Dehe Wang <dhwangwhu@whu.edu.cn>
* Miaomiao Wen <wenmiaomiao@whu.edu.cn>
* Kun Wang <wangk05@whu.edu.cn>
* Yu Zhou <yu.zhou@whu.edu.cn>
