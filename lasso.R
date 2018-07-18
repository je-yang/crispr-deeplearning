install.packages("glmnet")
library(glmnet)

seq_file = '/Users/JYang/Desktop/Stanford/modified_sequences.csv'
seq <- read.csv(seq_file, header = TRUE, stringsAsFactors = FALSE)
sapply(seq, mode)
seq_char <- as.character(seq$x)
seq_df_splits <- unlist(strsplit(seq_char, split = ""))
seq_matrix <- matrix(seq_df_splits , ncol = 38 , byrow = TRUE )
seq_df <- as.data.frame(seq_matrix)
