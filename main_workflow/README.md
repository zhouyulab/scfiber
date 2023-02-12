# Workflows

The analysis pipelines are implemented and managed using Snakemake.

- annotation.sm.rules: Generating gene annotations
- ATAC.sm.rules: ATAC-seq peak calling and motif analysis
- BuildNovelAnn.sm.rules: Generating novel gene annotations
- MetaInfo.sm.rules: Information of meta-data
- Nanopore.sm.rules: Data processing of SCAN-seq data
- RALF_RNA_seq.smk: Data processing of RALF RNA-seq data
- RNA_seq.sm.rule: Data processing of public Ovule RNA-seq 
- SingleCell3.sm.rules: UMAP, marker genes, DEGs and MEGs analysis for scRNA-seq data
- SingleCellPreTreatment3.sm.rules: Data processing and integration of scRNA-seq data

# Scripts

- ATAC_RNA.R: Association analysis between gene expression level and ATAC-seq signal level
- ATAC.profile.R: Signal profile of ATAC-seq data
- compute_pseudotime_enriched_gene.R: Finding pseudotime enriched genes
- compute_sc_tp_bin.R: Correlation analysis between scRNA-seq data and public ovule RNA-seq data
- kmer_graph_enrich.R: Identifying cell specific 8-mer motifs
- kmer_cluster2gene.R: Finding directional and non-directional cell specific 8-mer motifs
- plot_DE_gene_topGO.R: GO-TERM analysis of DEGs
- plot_kmer_enrich_motif.R: Figure plots of cell specific 8-mer motifs
- plot_marker_gene_num_topGO.R: GO-TERM analysis of marker genes
- plot_pseudotime_enriched_gene.R: Figure plots of expression pattern of pseudotime enriched genes
- SCAN_fiber_scRNA_cor.R: Correlation analysis among scRNA-seq, SCAN-seq, and Fiber RNA-seq
- scATAC_preprocess.R: UMAP and tSNE embedding analysis for scATAC-seq data
- scATAC_DE.R: Finding DARs and DAGs
