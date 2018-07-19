source("https://bioconductor.org/biocLite.R") 
biocLite("DESeq2")
library(DESeq2)

file = "/Users/JYang/Desktop/Stanford/20140519_ricintilingFinal_readcounts.csv"
data = read.csv(file, header = TRUE, stringsAsFactors=FALSE, row.names = "name")

#get columns of interest 
#dCas9 VP64
listnames1 <- colnames(data)[grepl('dCas9.VP64', colnames(data))]
index1 <- match(listnames1, names(data))
index1 <- sort(c(index1))
dCas9_CP64 <- data[ , index1] 

#scFV VP64
listnames2 <- colnames(data)[grepl('scFV.VP64', colnames(data))]
index2 <- match(listnames2, names(data))
index2 <- sort(c(index2))
scFV_VP64 <- data[ , index2] 

#create combined dataframe
dataset <- cbind(dCas9_CP64, scFV_VP64)

#split new dataset to separate ricin and cycled
col_names <- colnames(dataset)
condition <- c()
for(i in 1:length(col_names)){
  if(grepl('ricin', col_names[i])){
    condition <- c(condition, 'ricin')
  }
  else{
    condition <- c(condition, 'cycled')
  }
}
condition

coldata <- cbind(col_names, condition)
counts <- dataset
coldata <- as.matrix(coldata, row.names = 'col_names')

head(counts,2)
coldata

all(rownames(coldata) %in% colnames(counts))
all(rownames(coldata) == colnames(counts))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds

dds1 <- DESeq(dds)
res <- results(dds1)
res

res <- results(dds1, name="condition_ricin_vs_cycled")
res <- results(dds1, contrast=c("condition","ricin","cycled"))
res

log2fc <- res$log2FoldChange
log2fc
