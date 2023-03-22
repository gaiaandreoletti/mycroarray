library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(ggplot2)
# BiocManager::install("affy")
BiocManager::install("oligo")
library(oligo)
library(limma)
library(Biobase)
library(Biostrings)
library(genefilter)
# BiocManager::install("annotationTools")
library(annotationTools)
#BiocManager::install("pdInfoBuilder")
library(pdInfoBuilder)
BiocManager::install("pd.ht.hg.u133a")
library(pd.ht.hg.u133a) #https://support.bioconductor.org/p/55779/
# library(pd.u133aaofav2)

# tutorial: http://homer.ucsd.edu/homer/basicTutorial/affymetrix.html

setwd("~/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/CD110/GSE174060/")
# specify the path on your computer where the folder that contains the CEL-files is located
celpath = ("/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/CD110/GSE174060/GSE174060_RAW")
# # import CEL files containing raw probe-level data into an R AffyBatch object
# data = ReadAffy(celfile.path=celpath)
# import CEL files containing raw probe-level data into an R AffyBatch object
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list, pkgname="pd.ht.hg.u133a")
ph = data@phenoData

#normalize the data
eset <- rma(data)

# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="data.txt")

#Adding Gene Annotation to Normalized Expression Output
# Strategy is to create data frame objects and merge them together - put expression info into a data frame
my_frame <- data.frame(exprs(eset))
exprSet.nologs = exprs(eset)

# Put annotation information in a data frame.  To get specific fields, use packageNameSYMBOL, where the caps part names the type of data you're after
# To get a list of available annotation information, run the packagename with () at the end, i.e. mogene20sttranscriptcluster()
# Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
heatmap(exprSet.nologs, main = "Normalized ME matrix")


Annot <- read.table(file = "GPL4685_noParents.an.txt", header = TRUE)
head(Annot)
# Merge data frames together (like a database table join)
all <- merge(Annot, my_frame, by.x="ProbeName", by.y=0, all=T)
head(all)
# Write out to a file:
write.table(all,file="data.ann.txt",sep="\t")

samples <- read.table("sample.txt", header = T, row.names = 1)
samples$Cel_file_name <- rownames(samples)
dim(samples)

x.mas5 <- call.exprs(x,"mas5") # Calculates expression values with MAS 5.0 method which is required for the next step!

########################
# data_counts_cases <- read.table(file = "star_deseq_raw_counts_batch6.txt", header = TRUE, row.names = 1)
data_counts <- read.table(file = "/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/CD110/GSE174060/data.txt", header = TRUE)
# header.true <- function(df) {
#   names(df) <- as.character(unlist(df[1,]))
#   df[-1,]
# }
# df1 <- header.true(data_counts)
data_counts <- na.omit(data_counts)
# rownames(data_counts) <- make.unique(data_counts$GeneSymbols)
head(data_counts,2)
# data_counts$GeneSymbols <- NULL
dim(data_counts)

samples <- read.table("/Users/gandreoletti/Library/CloudStorage/OneDrive-GraphiteBio/CD110/GSE174060/sample.txt", header = T)
rownames(samples) <- samples$sample
dim(samples)
head(samples )
samples$condition <- as.factor(samples$condition)

library(DESeq2)
ncol(data_counts) == nrow(samples)
data_counts <- data_counts[,order(colnames(data_counts))]
samples <- samples[order(rownames(samples)),]
colnames(data_counts) == rownames(samples)


dds <- DESeqDataSetFromMatrix(countData = round(data_counts),
                              colData = samples,
                              design = ~ condition)


cts = counts(dds)
# geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
# dds = estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
dds<-estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
dds <- DESeq(dds, test = "Wald")
summary(dds)
vsd <- vst(dds, blind=FALSE, fitType='local')

d <- plotCounts(dds, gene="MPL", intgroup="Diagnosis", returnData=TRUE)

ggboxplot(d, x = "Diagnosis", y = "count", 
          color = "Diagnosis", palette =c("#00AFBB", "#E7B800", "#D55E00","#CC79A7","#F0E442", "#0072B2"),
          add = "jitter", shape = "Diagnosis") + ggtitle("CD110 expression across samples") + theme_bw() +
  stat_compare_means(label.y = 15)                   # Add global p-value
ggsave("CD110_expGSE54646_CEL.pdf", height=6, width=5.5)

d$Diagnosis <- factor(d$Diagnosis, levels=c("Normal", "ET", "MF","PV","JAK2_KO","HEL_Control"))
ggboxplot(d, x = "Diagnosis", y = "count", 
          color = "Diagnosis", palette =c("#00AFBB", "#E7B800", "#D55E00","#CC79A7","#F0E442", "#0072B2"),
          add = "jitter", shape = "Diagnosis") + ggtitle("CD110 expression across samples") + theme_bw() +
   stat_compare_means(comparisons = list(c("Normal", "ET"), c("Normal", "MF"), c("Normal", "PV"),c("Normal", "JAK2_KO")))
ggsave("CD110_expGSE54646_CEL_byGroup.pdf", height=6, width=7)

