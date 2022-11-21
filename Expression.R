library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

library("optparse")

option_list = list(
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="Write the metadata file name it must be .csv file", metavar="character"),
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="Path to the ballgown directory", metavar="character"),
  make_option(c("-c", "--condition"), type="character", default=NULL, 
              help="Condition column name from csv file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("All four argument must be supplied (input file).n", call.=FALSE)
}
# change this to the directory that contains all the StringTie results
# load the sample information
pheno_data <- read.csv(opt$metadata, header = TRUE)
# create a ballgown object
bg_data <- ballgown(dataDir = opt$path,
                    samplePattern = "SRR",
                    pData = pheno_data)
#rowVars calculated vatiance of each row
bg_filt = subset(bg_data,"rowVars(texpr(bg_data)) >1",genomesubset=TRUE) #subset subset only specific samples og genomic locations
bg_table = texpr(bg_filt, 'all') #texpr (extract transcrip-level expressions from bg objects )
bg_gene_names = unique(bg_table[, 9:10])  #takes unique values(all rows) from column 9 and 10 which gene id and gene name

gene_expression = as.data.frame(gexpr(bg_filt))
colnames(gene_expression) <- gsub(pattern = "FPKM.", replacement = "", x = names(gene_expression))
name <- as.matrix(geneNames(bg_filt))
dim(name)
#results_transcripts <- stattest(bg_filt, feature="transcript", covariate="phenotype", getFC=TRUE, meas="FPKM", adjustvars = NULL)
results_genes <- stattest(bg_filt, feature="gene", covariate=opt$condition, log = TRUE, meas="FPKM",  adjustvars = NULL, getFC = TRUE)

write.csv(bg_gene_names, "bg_gene_names.csv")
write.csv(bg_table, "bg_table.csv")
write.csv(results_genes, "results_genes.csv")
#results_transcriptsmer = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), transcriptNames=ballgown::transcriptNames(bg_filt), results_transcripts)
results_genes <- arrange(results_genes,pval)

#results_transcripts$logFC <- log2(results_transcripts$fc)
results_genes$logFC <- log2(results_genes$fc) #taking log of Fold Chanage Values

#write.csv(results_transcripts, "transcript_results.csv",row.names=FALSE)
write.csv(results_genes, "gene_results_log_fc_p_value_arranged.csv", row.names=FALSE) 

#subset(results_transcripts,results_transcripts$qval<0.05)
#subset(results_genes,results_genes$pval<0.05)

# subset of result transcripts where pvalue is less than 0.05
#tra <- subset(results_transcripts,results_transcripts$pval<0.05)
gen <- subset(results_genes,results_genes$pval<0.05) #confirm it first

#write.csv(tra, "filtered transcripts.csv")
write.csv(gen, "filtered_genes_pval_0.05.csv") #confirm it first

#####################Box plot#######################
tropical= c('darkorange', 'dodgerblue',
            'hotpink', 'limegreen', 'yellow')
palette(tropical)
#extract transcript-level expression measurements from ballgown objects
fpkm = texpr(bg_filt,meas="FPKM")
fpkm = log2(fpkm+1)
# 1. Open jpeg file
jpeg("boxplot.jpg", width = 350, height = 350)
boxplot(fpkm,col=as.numeric(pheno_data$phenotype),las=2,ylab='log2(FPKM+1)')
# 3. Close the file
dev.off()
write.csv(gene_expression, "raw_gene_expression_results.csv", row.names=FALSE)
write.csv(merge(results_genes, gene_expression, by=0, all=T), "merged_all_results.csv", row.names=TRUE)
#write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "final_logFCgene_results.csv", row.names=FALSE) 
transcript_gene_table = indexes(bg_filt)$t2g
head(transcript_gene_table)
#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
length(row.names(transcript_gene_table)) 
length(unique(transcript_gene_table[,"g_id"]))

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
jpeg("transcript_count.jpg", width = 350, height = 350)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)
dev.off()
#Plot #2 - the distribution of transcript sizes as a histogram
full_table <- texpr(bg_filt , 'all')
jpeg("transcript_length.jpg", width = 350, height = 350)
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")
dev.off()

#########################################
sig=which(results_genes$pval<0.05)
jpeg("DEGs_distribution.jpg", width = 350, height = 350)
hist(results_genes[sig,"logFC"], breaks=50, col="seagreen", xlab="log2(Fold change Disease vs. Normal)", main="Distribution of DEGs Values")
abline(v=-1, col="black", lwd=2, lty=2)
abline(v=1, col="black", lwd=2, lty=2)
legend("topright", "Fold-change > 4", lwd=2, lty=2)
dev.off()
library(EnhancedVolcano)

indices <- match(results_genes$id, texpr(bg_data, 'all')$gene_id)
gene_names_for_result <- texpr(bg_data, 'all')$gene_name[indices]
final_data_Frame <- data.frame(geneNames=gene_names_for_result, results_genes)
write.csv(final_data_Frame,"deg_result_genes_with_gene_names.csv")
gene_expression$GeneName <- gene_names_for_result
gene_expression <- subset(gene_expression, GeneName!=".")
write.csv(gene_expression, "FPKM_expression_NvP.csv", row.names = FALSE)
final_without_dots <- subset(final_data_Frame, geneNames!=".")
jpeg("Volcano_plot.jpg", width = 350, height = 350)
EnhancedVolcano(final_without_dots,
                lab = final_without_dots$geneNames,
                x = "logFC",
                y = "pval",
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()
