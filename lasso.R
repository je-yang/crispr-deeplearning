install.packages("glmnet")
library(glmnet)

#prepare data for lasso regression model 

#turn sequence position into categorical variables 
seq_file = '/Users/JYang/Desktop/Stanford/modified_sequences.csv'
seq <- read.csv(seq_file, header = TRUE, stringsAsFactors = FALSE)
sapply(seq, mode)
seq_char <- as.character(seq$x)
seq_df_splits <- unlist(strsplit(seq_char, split = ""))
seq_matrix <- matrix(seq_df_splits , ncol = 38 , byrow = TRUE )
seq_df <- as.data.frame(seq_matrix)

log_fold_changes_file = "/Users/JYang/Desktop/Stanford/log2fc.csv"
log_fold_changes_df <- read.csv(log_fold_changes_file, header = TRUE, stringsAsFactors = FALSE)


#bam_df from modify_sequences.R
#indices to remove
remove_idx <- which(bam_df$flag == 4)
seqs <- read.csv(file = "/Users/JYang/Desktop/Stanford/modified_sequences_unfiltered.csv", header = TRUE, stringsAsFactors = FALSE)
remove_idx2 <- which(strapply(seqs$sequence, "(..).....$", simplify = TRUE) != "GG")
#remove all unmapped
bam_df_mod <- bam_df[-remove_idx, ]
#remove ones without 'GG' in 32/33 position
bam_df_mod2 <- bam_df_mod[-remove_idx2, ]
#merge
df <- cbind(bam_df_mod2$qname, seq_df)
colnames(df)[1] <- "qname"
df_merge <- merge(df, log_fold_changes_df, by = 'qname')
final_df <- df_merge[complete.cases(df_merge), ]
final_df[40]
final_df <- final_df[, -c(40)]
#make all numeric columns absolute value
#final_df_abs <- cbind(final_df[, -c(40, 41, 42, 43, 44, 45)], abs(final_df[, c(40, 41, 42, 43, 44, 45)]))

#bin <- rbinom(nrow(final_df), 1, 0.5)

test_df <- final_df

###########


bassik_file_name = "/Users/JYang/Desktop/Stanford/bassik_table.csv"
bassik_file <- read.csv(bassik_file_name, header = TRUE, stringsAsFactors = FALSE)

names_list_file_name = "/Users/JYang/Desktop/Stanford/renamed_gene_counts.csv"
names_list <- read.csv(names_list_file_name, header = TRUE, stringsAsFactors = FALSE)

bin_val <- bassik_file$CRISPRa_phenotype
bin_val[bin_val < 0] <- -1 
bin_val[bin_val > 0] <- 1 
gene_bin <- cbind(bassik_file$Gene, bin_val)

library(stringr)
test_df_mod <- test_df
names <- test_df_mod$qname
#print(length(names))
#print(str_extract(names[1], "[^_]+"))
nameslist <- c()
for(i in 1:length(names)){
  name_gene = str_extract(names[i], "[^_]+")
  nameslist <- c(nameslist, name_gene)
}

test_df <- cbind(nameslist, test_df_mod)

colnames(test_df_mod)[1] <- "Gene"
colnames(gene_bin)[1] <- "Gene"

gene_bin_merge <- merge(gene_bin, test_df_mod, by = "Gene")
write.csv(gene_bin_merge, file = "/Users/JYang/Desktop/Stanford/gene_effects_data.csv")

data.frame(gene_bin_merge$bin_val, gene_bin_merge$log2fc)


###########

test_df <- gene_bin_merge

xfactors <- model.matrix(test_df$bin_val ~ test_df$V1 + test_df$V2 + test_df$V3 + test_df$V4 + test_df$V5 + test_df$V6
                         + test_df$V7 + test_df$V8 + test_df$V9 + test_df$V10 + test_df$V11 + test_df$V12
                         + test_df$V13 + test_df$V14 + test_df$V15 + test_df$V16 + test_df$V17 + test_df$V18
                         + test_df$V19 + test_df$V20 + test_df$V21 + test_df$V22 + test_df$V23 + test_df$V24
                         + test_df$V25 + test_df$V26 + test_df$V27 + test_df$V28 + test_df$V29 + test_df$V30
                         + test_df$V31 + test_df$V34 + test_df$V35 + test_df$V36 + test_df$V37 + test_df$V38)[, -1]
x        <- as.matrix(data.frame(test_df$log2fc, xfactors))

# Note alpha=1 for lasso only and can blend with ridge penalty down to
# alpha=0 ridge only.
#test_mat <- as.matrix(test_df[,-c(1)])
glmmod <- glmnet(x, y=test_df$bin_val, alpha=1, family="binomial")

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")
glmmod
coef(glmmod)[, 5]
cv.glmmod <- cv.glmnet(x, y=factor(test_df$bin_val), alpha=1, family="binomial")
plot(cv.glmmod)
(best.lambda <- cv.glmmod$lambda.min)

test_df_1 <- data.frame(test_df$Gene, test_df$log2fc)
keys <- colnames(test_df_1)[!grepl('test_df.log2fc',colnames(test_df_1))]
library(data.table)
X <- as.data.table(test_df_1)
X[,list(mm= mean(test_df.log2fc)),keys]
list <- X[,list(mm= mean(test_df.log2fc)),keys]
colnames(list)[1] <- "Gene"
list
check <- merge(list, test_df, by = "Gene")
check.original <-check
check$mm[check$mm < 0] <- 0 
check$mm[check$mm > 0] <- 1 
View(check)
sum(check$mm != check$bin_val)


list.med <- X[,list(mm= median(test_df.log2fc)),keys]
colnames(list)[1] <- "Gene"
list.med
check <- merge(list, test_df, by = "Gene")
check.original <-check
check$mm[check$mm < 0] <- 0 
check$mm[check$mm > 0] <- 1 
View(check)
sum(check$mm != check$bin_val)
