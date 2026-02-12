#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(optparse))
suppressMessages(library(clustree))

###Script to plot clustree graph from clustree package
option_list <- list(
  make_option(c("-o", "--output_dir"), action = "store", type = "character",
              default = getwd(), help = "Directory where the output files will be saved.[default: current working directory]"),
  make_option(c("-v", "--version"),
              action = "store_true", default = FALSE,
              help = "Print the list of loaded package versions and exit.")
)

parser <- OptionParser(usage = "%prog [options] csv_file",
                       option_list = option_list, prog = "clustree",
                       description = "Visualize clustree plot using the clustree R package.'csv_file' is the path of the clustering_ID for each cell at each resolution.")

arguments <- parse_args(parser, args <- commandArgs(trailingOnly = TRUE),
                        positional_arguments = TRUE)
opt <- arguments$options
version <- opt$version
output_dir <- opt$output_dir
csv_file <- arguments$args[1]

if(version){
  # Printing package versions
  x = sessionInfo()
  cat(paste(" "," "," R: ", x$R.version$major,".", x$R.version$minor, "\n", sep=""))
  for(i in 1:length(x$otherPkgs)){
    cat(paste(" "," "," ", x$otherPkgs[[i]]$Package,": ", x$otherPkgs[[i]]$Version, "\n", sep=""))
  }
  quit()
}

# Importing csv file
df <- read.csv(csv_file, h = TRUE, row.names = 1, check.names = FALSE)
df_filtered <- df[, startsWith(names(df), "leiden")]


# Defining the output file path in the specified directory
output_file <- paste(output_dir, "/clustree_plot.pdf", sep = "")

# Check if more than one single resolution has been computed
if(!is.vector(df_filtered)){
  #Save clustree plot
  pdf(output_file, width=10, height=12)
  print(clustree(df_filtered, prefix = "leiden_", suffix = ""))
  dev.off()
} else {
  #Save warning message
  pdf(output_file, width=10, height=12)
  plot.new()
  par(mar = c(0, 0, 1.1, 0))
  text(x=0.5, y=0.5, labels="One single resolution was found in the considered object. \nClustree not computed.", cex=1.5)
  dev.off()
}
