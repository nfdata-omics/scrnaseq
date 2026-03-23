#!/usr/bin/env Rscript
options(warn=-1)
suppressMessages(library("optparse"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("ggplot2"))

### Script to fit a DGE model using DESeq2
option_list <- list(
  make_option(c("-F", "--FDR"), action="store", type="double", default=0.05, help="The FDR cutoff to define significant genes."),
  make_option(c("-s","--suffix"), action="store", type="character", default="", help="Suffix to append to the output filenames, e.g. deseq2_toptable.variable_num_vs_denom_SUFFIX.txt."),
  make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print the list of loaded package versions and exit.")
)
### change logFC names
parser<-OptionParser(usage = "%prog [options] input_model contrast",
                     option_list = option_list, prog = "dge_deseq2_results",
                     description = "
                     Extract results from a DGE model fitted with DESeq2.
                     'input_model' is the path of the DESeq2 object containing the already fitted model.
                     'contrast' is a string of format variable:test:reference, with words separated by colon (:), where
                                - 'variable' is the name of the variable (metadata column) on which the comparison will be performed (e.g. treatment). It must be one of the variables previously included in the model formula;
                                - 'test' is the name of the factor level to be tested (the numerator of the log2 Fold Change, e.g. treat). If multiple levels should be aggregated, they must be separated by '/' with no blank spaces (e.g. treat1/treat2/treat3);
                                - 'reference' is the name of the factor level to be used as reference (the denominator of the log2 Fold Change, e.g. ctrl). If multiple levels should be aggregated, they must be separated by '/' with no blank spaces (e.g. ctrl1/ctrl2)."
)

arguments <- parse_args(parser, args <- commandArgs(trailingOnly=TRUE), positional_arguments = TRUE)
opt <- arguments$options
version = opt$version

if (length(arguments$args)!=2 & !version) {
  stop("Two arguments must be supplied (input_model and contrast)", call.=FALSE)
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

input_model = arguments$args[1]
contrast = str_split_1(arguments$args[2], ":")
variable = contrast[1]
num = str_split_1(contrast[2], "/")
n_num = length(num)
numNames = gsub("/", "_", contrast[2])
denom = str_split_1(contrast[3], "/")
n_denom = length(denom)
denomNames = gsub("/", "_", contrast[3])
fdr = opt$FDR
suffix = opt$suffix

# Importing dge object
dds = get(load(input_model))

# Extracting results
if(n_num>1 | n_denom>1){

  # Custom list A1,A2,A3... vs B1,B2... contrast
  num_list = paste(variable, num, sep="")
  num_list = num_list[num_list%in%resultsNames(dds)]
  n_num = length(num_list)
  numNames = gsub(variable, "", paste(num_list, collapse="_"))
  denom_list = paste(variable, denom, sep="")
  denom_list = denom_list[denom_list%in%resultsNames(dds)]
  n_denom = length(denom_list)
  denomNames = gsub(variable, "", paste(denom_list, collapse="_"))
  my_res = results(dds, contrast=list(num_list, denom_list), listValues=c(1/n_num, -1/n_denom), independentFiltering=TRUE, cooksCutoff=FALSE, alpha=fdr)
  # I turn off cooksCutoff for outlier detection, but in the model fitting function there was minRepforReplace=7

} else {

  # Direct A vs B contrast
  my_res = results(dds, contrast=c(variable, num, denom), independentFiltering=TRUE, cooksCutoff=FALSE, alpha=fdr)
  # I turn off cooksCutoff for outlier detection, but in the model fitting function there was minRepforReplace=7

}

capture.output(summary(my_res, alpha=fdr), file=paste("deseq2_summary.",variable,"_",numNames,"_vs_",denomNames,suffix,".txt", sep=""))


pdf(paste("MAplot.",variable,"_",numNames,"_vs_",denomNames,".pdf", sep=""))
DESeq2::plotMA(my_res, ylim=c(-6,6), alpha=fdr, main=paste(numNames," vs ",denomNames, sep=""))
dev.off()


# Saving toptable
n_info = (which(colnames(rowData(dds))=="baseMean")-1)
toptable = cbind(my_res[,c("baseMean","log2FoldChange","pvalue","padj")], rep("not_sign", nrow(my_res)),rowData(dds)[,1:n_info])
colnames(toptable)[5:ncol(toptable)] = c("Significance", colnames(rowData(dds))[1:n_info])
toptable$Significance[!is.na(toptable$padj) & toptable$padj < fdr & toptable$log2FoldChange > 0] = "up"
toptable$Significance[!is.na(toptable$padj) & toptable$padj < fdr & toptable$log2FoldChange <= 0] = "down"
toptable = toptable[order(toptable$pvalue, decreasing=F),]
# Excluding the case in which rowData(rna_exp) ARE exactly the rownames
if(sum(toptable[,6]==rownames(toptable)) == nrow(toptable)) {toptable = toptable[,-6]}
toptable = cbind(rownames(toptable), toptable)
colnames(toptable)[1] = "gene_name"
write.table(toptable, paste("deseq2_toptable.",variable,"_",numNames,"_vs_",denomNames,suffix,".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")



# Histogram of nominal pvalues
p = ggplot(toptable, aes(x=pvalue)) + geom_histogram(binwidth=0.05, color="blue", fill="skyblue1") + theme_classic()
p = p + labs(title=paste("Nominal pvalue distribution - ", variable, ": ", numNames," vs ",denomNames, sep=""), x="nominal pvalue", y="frequency") + theme(title=element_text(size=13), axis.text=element_text(size=13))
pdf(paste("pvalue_hist.",variable,"_",numNames,"_vs_",denomNames,suffix,".pdf",sep=""), width=8, height=8)
print(p)
dev.off()


# Volcano plot
colvol = rep("gray", nrow(toptable))
colvol[toptable$padj<fdr & toptable$log2FoldChange>0] = "red2"
colvol[toptable$padj<fdr & toptable$log2FoldChange<0] = "royalblue"
logpv = -log10(toptable$pvalue)
new_data = cbind(toptable[,c(3:5)], logpv, colvol)
colnames(new_data)=c("logFC", "pvalue", "FDR", "logpv", "sign")
subt = paste("Significance level was defined using the cutoff FDR<",fdr, sep="")
new_data = new_data[!is.na(new_data$FDR),]

pval_cutoff = tail(new_data$pvalue[new_data$FDR<fdr], 1) # the raw pvalue corresponding to the last significant FDR
p = ggplot(new_data, aes(x=logFC, y=logpv, color=sign)) + geom_point(shape=19, show.legend=F)
p = p + scale_color_manual(values = c("gray" = "gray", "red2" = "red2", "royalblue" = "royalblue"))
p = p + geom_hline(yintercept=-log10(pval_cutoff), linetype="dashed", col=1) + geom_vline(xintercept=c(-1,1), linetype="dashed", col=1)
p = p + labs(title=paste("Volcano plot - ",variable, ": ", numNames," vs ",denomNames, sep=""), subtitle=subt, x="log2FC", y="- log10 (pvalue)") + theme_classic()
p = p + theme(title=element_text(size=20), axis.text=element_text(size=15))
pdf(paste("volcano.",variable,"_",numNames,"_vs_",denomNames,suffix, ".pdf", sep=""), width=10, height=8)
print(p)
dev.off()

########################

w <- warnings()
if (!is.null(w)) {
  writeLines(capture.output(print(w)), con = stderr())
}
