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
remove_idx2 <- which(strapply(seqs, "(..).....$", simplify = TRUE) != "GG")
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
