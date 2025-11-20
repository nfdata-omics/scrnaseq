#!/usr/bin/env Rscript

options(warn = -1)
suppressMessages(library(optparse))
suppressMessages(library(scDblFinder))
suppressMessages(library(SingleCellExperiment))

###Script to compute doublets from scdblfinder package
option_list <- list(make_option(c("-o", "--output_dir"), action = "store", type = "character", default = getwd(), help = "Directory where the output files will be saved. [default: current working directory]"),
            make_option(c("-v", "--version"), action = "store_true", default = FALSE, help = "Print the list of loaded package versions and exit.")
            )
parser <- OptionParser(usage = "%prog [options] *.sce", option_list = option_list, prog = "scDblFinder", description = "Compute doublets score for each cell in each sample using the scDblFinder R package.'*.sce' is the path of the concatenated count matrix as SingleCellExperiment object.")
arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
opt <- arguments$options
version <- opt$version
output_dir <- opt$output_dir

if(version){
  # Printing package versions
  x = sessionInfo()
  cat(paste(" "," "," R: ", x$R.version$major,".", x$R.version$minor, "\n", sep=""))
  for(i in 1:length(x$otherPkgs)){
    cat(paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, "\n", sep=""))
  }
  quit()
}

sce_object <- arguments$args[1]

# Importing SingleCellExperiment object
sce <- readRDS(sce_object)
print(sce)

# Remove columns starting with "fastq"
fastq_cols <- grep("^fastq", colnames(colData(sce)), value = TRUE)
colData(sce) <- colData(sce)[, !(colnames(colData(sce)) %in% fastq_cols)]

names(assays(sce)) <- "counts"

#Compute doublets
sce <- scDblFinder(sce, sample = "sample")

# Defining the output file path in the specified directory
cell_annotation <- as.data.frame(sce@colData)
cell_annotation = cbind(rownames(cell_annotation), cell_annotation)

#Save csv file
output_file <- paste(output_dir, "/doublets_score.csv", sep = "")
write.csv(cell_annotation, output_file, quote=F, row.names=F, col.names=TRUE)
