#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("openxlsx"))


### Script to perform over-representation enrichment on marker genes using clusterProfiler
option_list <- list(
  make_option(c("-p", "--pval_adj"), action="store", type="double", default=0.05, help="The pvalue-adjusted cutoff to define significant genes. [default \"%default\"]"),
  make_option(c("-l", "--logfc"), action="store", type="double", default=0.25, help="The log2FC cutoff to define significant genes (in combination with adjusted pvalue). [default \"%default\"]"),
  make_option(c("-b", "--background"), action="store_true", default=FALSE, help="Whether to consider the full list of marker genes from all the clusters as the experiment background"),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser<-OptionParser(usage = "%prog [options] marker_genes resolution gmt_file",
                     option_list = option_list, prog = "markers_clusterProfiler",
                     description = "
                     Perform over-representation analysis using the clusterProfiler R package. 
                     'marker_genes' is the path of the .xlsx table from which marker genes are extracted. It must be produced by the clustering step of the scrnaseq pipeline.
                     'resolution' refers to the selected Leiden resolution. The corresponding clustering must have been already calculated.
                     'gmt_file' is the path of a .gmt file containing pathways to evaluate. Each pathway must be reported in a different row, with genes separated by tabs.")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=3 & !version) {
  stop("Three arguments must be supplied (marker_genes, resolution and gmt_file)", call.=FALSE)
}

if(version){
  # Printing package versions
  x = sessionInfo()
  cat(paste(" "," "," R: ", x$R.version$major,".", x$R.version$minor, "\n", sep=""))
  for(i in 1:length(x$otherPkgs)){
    cat(paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, "\n", sep=""))
  }
  quit()
}


marker_genes = arguments$args[1]
res = arguments$args[2]
gmt_file = arguments$args[3]
pval_adj = opt$pval_adj
logfc = opt$logfc
use_background = opt$background


# Import gmt file
my_gmt = read.gmt(gmt_file)
# Uppercase to guarantee consistency also from Enrichr-downloaded gmt
my_gmt$gene = toupper(my_gmt$gene)

# Importing marker genes for the considered resolution
my_table = read.xlsx(marker_genes, sheet=paste0("Leiden_",res))
clusters = names(table(my_table$group))
background_genes = toupper(unique(my_table$names))

enrich_list = list()
for(i in clusters){
  my_genes = my_table$names[my_table$group==i & my_table$pvals_adj<pval_adj & my_table$logfoldchanges>logfc]
  my_genes = toupper(my_genes)
  if(length(my_genes)>5){
    if(use_background){
      my_enrich = tryCatch({
        enricher(gene=my_genes, universe=background_genes, maxGSSize=1000, TERM2GENE=my_gmt)},
        # BgRatio = M/N, with M = #genes in the selected pathway that are also in the background; with N = #unique genes in all the gmt (all pathways) that are also in the background
        # GeneRatio = M/N, with M = #my input genes that are also in the selected pathway; with N = #my input genes that are also in all the gmt (all pathways)
        error = function(e) {message(e$message)
        })
    } else {
      my_enrich = tryCatch({
        enricher(gene=my_genes, maxGSSize=1000, TERM2GENE=my_gmt)},
        error = function(e) {message(e$message)
        })
    }
    if(!is.null(my_enrich)) {
      enrich_list[[paste0("cl_",i)]] = my_enrich@result
    } else {
      no_enrichment = data.frame(empty="No enrichment was found.")
      enrich_list[[paste0("cl_",i)]] = no_enrichment
    }
  } else {
    no_enrichment = data.frame(empty=paste("Enrichment analysis was not performed. Less than 5 marker genes were retrieved using the imposed significance cutoff: pval_adj<", pval_adj, " and logFC>", logfc, sep=""))
    enrich_list[[paste0("cl_",i)]] = no_enrichment
  }
}

collection = substring(basename(gmt_file), 1, nchar(basename(gmt_file))-4)
outname = paste("enrich_Leiden_", res, "_", collection, ".xlsx", sep="")
write.xlsx(enrich_list, file=outname)





########################

w=warnings()
sink(stderr())
if(!is.null(w)){
  print(w)
}
sink()
