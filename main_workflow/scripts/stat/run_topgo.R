library(topGO)
library(readr)
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--gene-list", required=TRUE, type="character", dest = "gene_list", metavar="gene_list.txt")
parser$add_argument("--background-list", required=TRUE, type="character", dest = "bg_list", metavar="bg_list.txt")
parser$add_argument("--GO-db", required=TRUE, type="character", dest = "GO_db", metavar="GO.map")
parser$add_argument("--top-nodes", required=FALSE, type="integer", dest = "top_nodes", metavar="top_nodes", default=10)
parser$add_argument("--output", required=TRUE, type="character", dest = "output", metavar="output")
args <- commandArgs(TRUE)
args <- parser$parse_args(args) 

gene_list <- as.character(unlist(read_csv(args$gene_list, col_names = FALSE)))
bg_list <- as.character(unlist(read_csv(args$bg_list, col_names = FALSE)))
geneID2GO <- readMappings(file = args$GO_db)

geneList <- factor(as.integer(bg_list %in% gene_list), levels=c(0, 1))
names(geneList) <- bg_list
print(sprintf("Gene list: %d, Background gene list: %d", length(gene_list), length(geneList)))

MF_GOdata <- new("topGOdata", ontology="MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
BP_GOdata <- new("topGOdata", ontology="BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
CC_GOdata <- new("topGOdata", ontology="CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

go_data2result <- function(go_data, f_out, ontology, geneList, top_nodes=10){
  resultFisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(go_data, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(go_data, algorithm = "elim", statistic = "ks")
  allRes <- GenTable(go_data, classicFisher = resultFisher, 
                     classicKS = resultKS, elimKS = resultKS.elim, 
                     orderBy = "classicFisher", ranksOf = "classicFisher", 
                     topNodes = 100)
  allRes <- allRes[allRes$Significant > 0,]
  res <- cbind(data.frame(Ontology=ontology), allRes)
  res$Gids <- sapply(res$GO.ID, function(x){return(paste(intersect(geneList, genesInTerm(go_data, x)[[1]]),collapse = " "))})
#  printGraph(go_data, resultFisher, firstSigNodes = top_nodes, useInfo = "all", fn.prefix = file.path(f_out, ontology))
  return(res)
}

MF_ref <- go_data2result(MF_GOdata, args$output, "MF", gene_list, args$top_nodes)
BP_ref <- go_data2result(BP_GOdata, args$output, "BP", gene_list, args$top_nodes)
CC_ref <- go_data2result(CC_GOdata, args$output, "CC", gene_list, args$top_nodes)
df <- rbind(MF_ref, BP_ref, CC_ref)
write_tsv(df, file.path(args$output, "topGO.tsv"))
