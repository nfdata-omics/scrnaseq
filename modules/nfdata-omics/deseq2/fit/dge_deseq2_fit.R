#!/usr/bin/env Rscript

options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("SummarizedExperiment"))

### Script to fit a DGE model using DESeq2
option_list <- list(
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output paths, e.g. deseq2_obj_SUFFIX.rds."),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)

parser <- OptionParser(usage = "%prog [options] counts metadata model_formula",
                     option_list = option_list, prog = "dge_deseq2_fit",
                     description = "
                     Fit a differential gene expression (DGE) model using DESeq2.
                     'counts' is the path of a .tsv file containing the expression data. This file should have been previously filtered to remove lowly expressed features and normalized for library size (e.g. using CPM normalization). The first column of the file should contain the feature IDs and the first row should contain the sample names.
                     'metadata' is the path of a .tsv file containing the samples metadata. The first column of the file should contain the sample names and the first row should contain the variable names.
                     'model_formula' is the formula used for the design of the DESeq2 model. It must start with a '~' and include all the biological and technical variables that should be accounted for (e.g. ~0+genotype+treatment+batch) with no blank spaces.
                                     If you plan to perform a comparison between aggregated combinations of levels (e.g. treat1/treat2/treat3 vs ctrl), the model formula must start with zero '~0+... ', otherwise zero could be omitted."
)

# Parse arguments
arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=3 & !version) {
  stop("Three arguments must be supplied (counts, metadata and model_formula)", call.=FALSE)
}

# Print package versions
if(version){
  x = sessionInfo()
  cat(paste(" "," "," R: ", x$R.version$major,".", x$R.version$minor, "\n", sep=""))
  for(i in 1:length(x$otherPkgs)){
    cat(paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, "\n", sep=""))
  }
  quit()
}

# Create variables with arguments
counts_path = arguments$args[1]
metadata_path = arguments$args[2]
model_formula = as.formula(arguments$args[3])
suffix = opt$suffix

# Import counts and set rownames
counts = read.delim(counts_path,
                    header = TRUE)
rownames(counts) = counts[,1]
counts = counts[,-1]

# Check if counts can be converted to numeric matrix. If not, stop the script and print an error message
if(!all(sapply(counts, is.numeric))){
  stop("Counts cannot be converted to numeric matrix.
       Please check the input file and make sure that are no gene names or other non-numeric values in the counts matrix.", call.=FALSE)
}

# Convert counts to numeric matrix
counts = as.matrix(counts)
mode(counts) = "numeric"

# Import metadata and set rownames
metadata = read.delim(metadata_path,
                      header = TRUE)
rownames(metadata) = metadata[,1]
metadata = metadata[,-1]

# Order metadata and counts by the same sample order
metadata = metadata[order(rownames(metadata)),]
counts = counts[,order(colnames(counts))]

# Define the DESeq2 object
dds = DESeqDataSetFromMatrix(countData=round(as.matrix(counts)),
                             colData=metadata,
                             rowData=counts,
                             design=model_formula)

colnames(rowData(dds)) = colnames(counts)

# DispEst plot for the whole dataset
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)

# Save the plot as a PDF file
pdf(paste("DispEst_plot",suffix,".pdf", sep=""), width=8, height=8)
plotDispEsts(dds)
dev.off()

# Fit models
dds = DESeq(object=dds,
            test="Wald",
            fitType="parametric",
            betaPrior=FALSE,
            minReplicatesForReplace=7,
            parallel=F)

save(dds, file=paste("deseq2_obj", suffix, ".rds", sep=""))

########################

w <- warnings()
if (!is.null(w)) {
  writeLines(capture.output(print(w)), con = stderr())
}
