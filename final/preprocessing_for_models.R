install.packages("glmnet")
library(glmnet)

####################################
#FILES:
#
#Input:
#modified_sequences.csv (from modify_sequences.R)
#log2fc.csv (from log_fold_chage_DESeq.R)
#modified_sequences_unfiltered.csv
#bassik_table.csv
#renamed_gene_counts.csv
#GilbertTreatedVsUntreatedGeneScores.csv
#gene_effect_means.csv
#
#Output:
#lasso_data.csv

####################################
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


####################################################
####################################################
####################################################

#Check processed data
#get density plot of fold change distribution

gene_scores_file = "/Users/JYang/Desktop/Stanford/GilbertTreatedVsUntreatedGeneScores.csv"
gene_scores_df <- read.csv(gene_scores_file, header = TRUE, stringsAsFactors = FALSE)

gene_effect_means_file = "/Users/JYang/Desktop/Stanford/gene_effect_means.csv"
gene_effect_means_df <- read.csv(gene_effect_means_file, header = TRUE, stringsAsFactors = FALSE)

gene_list <- gene_scores_df[gene_scores_df$Gene %in% gene_effect_means_df$Gene, ]
View(gene_list)

gene_list_combined <- merge(gene_list, gene_effect_means_df, by = "Gene")
summary(gene_list$FDR)
#mean is 0.2986
threshold = 0.5
genes_below_threshold <- gene_list_combined[gene_list_combined$FDR < threshold, ]
genes_above_threshold <- gene_list_combined[gene_list_combined$FDR >= threshold, ]
plot(density(genes_below_threshold$mm))
lines(density(genes_above_threshold$mm))


library(ggplot2)

#Sample data
dat <- data.frame(dens = c(genes_below_threshold$mm, genes_above_threshold$mm)
                  , lines = rep(c("<0.5", ">= 0.5")))

#Plot of log2fc density (distribution)
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5) + scale_x_continuous(limits = c(-0.5, 0.4)) +
  theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank(), plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) +
  ggtitle("Fold Change Distribution") + xlab("Average Log2 Fold Change") + ylab("Density") + labs(fill='Threshold')

#get final dataset for models

gene_list_mod <- gene_list
pos_vs_neg <- c()
for(i in 1:length(gene_list_mod$Gene)){
  if(gene_list_mod$positivePosteriorProb[i] > gene_list_mod$negativePosteriorProb[i]){
    pos_vs_neg <- c(pos_vs_neg, 1)
    print("neg")
  }
  if(gene_list_mod$positivePosteriorProb[i] < gene_list_mod$negativePosteriorProb[i]){
    pos_vs_neg <- c(pos_vs_neg, -1)
    print("pos")
  }
}
head(pos_vs_neg)
combined_check <- cbind(gene_list, pos_vs_neg)
gene_data_merge <- merge(gene_effect_means_df, combined_check, by = "Gene")
final_data_set <- merge(gene_data_merge, gene_bin_merge, by = "Gene")
final_data_set <- within(final_data_set, rm(log2fc))
colnames(final_data_set)[2] <- "avg_log2fc"
write.csv(final_data_set, row.names = FALSE, file= "/Users/JYang/Desktop/Stanford/lasso_data.csv")
final_data_set <- read.csv(file= "/Users/JYang/Desktop/Stanford/lasso_data.csv", header = TRUE)

