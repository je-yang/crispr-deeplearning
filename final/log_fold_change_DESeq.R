source("https://bioconductor.org/biocLite.R") 
biocLite("DESeq2")
library(DESeq2)

############################
#computing log2 fold changes with DESeq2

file = "/Users/JYang/Desktop/Stanford/20140519_ricintilingFinal_readcounts.csv"
data = read.csv(file, header = TRUE, stringsAsFactors=FALSE, row.names = "name")

#get columns of interesultst 
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

#check that row and columns match
all(rownames(coldata) %in% colnames(counts))
all(rownames(coldata) == colnames(counts))

library("DESeq2")
DESeq_data <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
DESeq_data

DESeq_data1 <- DESeq(DESeq_data)
results <- result(DESeq_data1)
results

results <- result(DESeq_data1, name="condition_ricin_vs_cycled")
results <- result(DESeq_data1, contrast=c("condition","ricin","cycled"))
results

log2fc <- cbind(rownames(results), results$log2FoldChange)
colnames(log2fc) <- c('qname', "log2fc")
head(log2fc)

log2fc_df <- data.frame(log2fc)

write.csv(log2fc, file = "/Users/JYang/Desktop/Stanford/log2fc.csv")
log2fc_df <- read.csv(file = "/Users/JYang/Desktop/Stanford/log2fc.csv", header =TRUE)

#check where mean and median are distributed
mean(log2fc_df$log2fc, na.rm=TRUE)
median(log2fc_df$log2fc, na.rm=TRUE)
