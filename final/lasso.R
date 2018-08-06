library(e1071)
library(rpart)
library(randomForest)
#install.packages('caret', dependencies = TRUE)
#install.packages("fields")
library(caret)
library(fields)
library(glmnet)

####################################
#FILES:
#
#Input:
#lasso_data.csv

############################

# load the CSV file from the local directory
filename <- "/Users/JYang/Desktop/Stanford/lasso_data.csv"
dataset_orig <- read.csv(filename, header=TRUE)

#randomize order of dataset
dataset <- dataset_orig[sample(nrow(dataset_orig)),]

#Split the loaded dataset into two
#80% of which we will use to train our models and 20% that we will hold back as a validation dataset
validation_index <- createDataPartition(dataset$avg_log2fc, p=0.80, list=FALSE)
# select 20% of the data for validation
validation <- dataset[-validation_index,]
# use the remaining 80% of data to training and testing the models
dataset <- dataset[validation_index,]
test_df <- dataset


x <- dummyVars(pos_vs_neg*avg_log2fc ~ Gene + V1 + V2 + V3 + V4 + V5 + V6
               + V7 + V8 + V9 + V10 + V11 + V12
               + V13 + V14 + V15 + V16 + V17 + V18
               + V19 + V20 + V21 + V22 + V23 + V24
               + V25 + V26 + V27 + V28 + V29 + V30
               + V31 + V34 + V35 + V36 + V37 + V38, data = test_df)

x = predict(x, newdata = test_df)

pred_test = test_df$pos_vs_neg*test_df$avg_log2fc

#validation set
pred_val <- as.numeric(validation$avg_log2fc)*as.numeric(validation$pos_vs_neg)
validation <- within(validation, rm(avg_log2fc, pos_vs_neg))
validation <- cbind(validation, pred_val)
validation <- within(validation, rm(positivePosteriorProb, negativePosteriorProb, FDR, qname, V32, V33))

newX <- dummyVars( ~ Gene + V1 + V2 + V3 + V4 + V5 + V6
                   + V7 + V8 + V9 + V10 + V11 + V12
                   + V13 + V14 + V15 + V16 + V17 + V18
                   + V19 + V20 + V21 + V22 + V23 + V24
                   + V25 + V26 + V27 + V28 + V29 + V30
                   + V31 + V34 + V35 + V36 + V37 + V38 , data = validation)

newX = predict(newX, newdata = validation)

glmmod <- glmnet(x, y=pred_test, alpha=1, intercept = FALSE)

# Plot variable coefficients vs. shrinkage parameter lambda.
plot(glmmod, xvar="lambda")
glmmod
cv.glmmod <- cv.glmnet(x, y=pred_test, alpha=1, nfolds=5)
plot(cv.glmmod)
best.lambda <- cv.glmmod$lambda.min
coef_table <- coef(cv.glmmod, s = "lambda.min")
coef_table_orig <- coef_table

#########Plot coef agaisnt phenotype values
gene_names_mm <- read.csv("/Users/JYang/Desktop/Stanford/gene_effect_means.csv", header = TRUE, stringsAsFactors = FALSE)
gene_names <- gene_names_mm[,1]
gene_means <- gene_names_mm[,2]
coef_list <- as.vector(coef_table_orig[-1])
coef_list <- head(coef_list, 48)
coef_list <- cbind(gene_names_mm, coef_list)
bassik_table <- read.csv("/Users/JYang/Desktop/Stanford/bassik_table.csv", header=TRUE)
coef_phen_list <- merge(bassik_table, coef_list, by="gene_names")
coef_phen_list$coef_list <- as.numeric(as.character(coef_phen_list$coef_list))


p <- ggplot(coef_phen_list, aes(CRISPRa_phenotype, coef_list)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.4, col = "black") +
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("CRISPRa Phenotype vs LASSO Coefficient") + xlab("CRISPRa Phenotype") + ylab("LASSO Coefficient")+ 
  ylim(-0.36,0.36) + xlim(-0.36,0.36)
p


g <- ggplot(coef_phen_list, aes(coef_list, mm, color=CRISPRa_phenotype)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.4) + 
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("LASSO Coefficient vs Log2FC") + xlab("LASSO Coefficient") + ylab("Log2FC")+ 
  ylim(-0.5,0.5) + xlim(-0.5,0.5)
g

filename <- "/Users/JYang/Desktop/Stanford/lasso_data.csv"
dataset_orig <- read.csv(filename, header=TRUE)
df_log2fc <- dataset_orig[!duplicated(dataset_orig$Gene), ]
coef_phen_log_list <- cbind(coef_phen_list, df_log2fc$avg_log2fc)
colnames(coef_phen_log_list)[6] <- "log2fc"

h <- ggplot(coef_phen_log_list, aes(coef_list, log2fc, color=CRISPRa_phenotype)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.4) + 
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("LASSO Coefficient vs Gene Mean") + xlab("LASSO Coefficient") + ylab("Gene Mean")+ 
  ylim(-0.5,0.5) + xlim(-0.5,0.5)
h


###############

z = matrix(as.vector(coef_table)[-1], nrow= 4)
base<- rownames(z)
pos <- colnames(z)
z<- data.frame(z)

#Test on validation set
#get mean squared error
fit_test<-predict(cv.glmmod, newx=newX, s="lambda.min")
plot(validation$pred_val, fit_test, col = "blue", pch=4)
mse.min <- min(cv.glmmod$cvm)
error <- validation$pred_val - fit_test
mse <- mean(error^2)

#install.packages("ggpubr")
library(ggpubr)

#get pearson correlation coefficient, comparing actual vs predicted
pearson_corr = cor(validation$pred_val, fit_test, method = "pearson") #"spearman"))

#plot actual vs predicted
library(ggplot2)
df <- data.frame(cbind(validation$pred_val, fit_test))
colnames(df)<-c('Actual', 'Predicted')
df$pc <- predict(prcomp(~validation$pred_val+fit_test, df))[,1]
p <- ggplot(df, aes(Actual, Predicted, color = pc)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.1) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") + 
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("LASSO Regression - Log Fold Change") + xlab("Actual") + ylab("Predicted") +
  annotate("text", x = c(-0.3), y=c(0.1), label = paste0("Pearson Correlation: ", "\n", pearson_corr))
p


################
#Get heat map of base, position and LASSO Coefficients

mat_mod <- data.frame(coef.name = dimnames(coef(glmmod))[[1]], coef.value = as.matrix(coef(glmmod)))

positions<- as.numeric(gsub("[^0-9\\.]", "", mat_mod$coef.name))
positions
mat_mod1 <- cbind(positions, mat_mod)
require(stringi)
mat_mod1$coef.name <- stri_sub(mat_mod1$coef.name, from=-1, to=-1)
mat_mod1 <- mat_mod1[c(-1),]

mat_mod1_wide <- reshape(mat_mod1, idvar = "positions", timevar = "coef.name", direction = "wide")

z = matrix(as.vector(coef_table)[-1], nrow= 4)
base<- rownames(z)
pos <- colnames(z)
z<- data.frame(z)

require(reshape2)
z_melted <- cbind(data.frame(matrix(as.vector(coef_table)[-1])),mat_mod1$positions,mat_mod1$coef.name)
colnames(z_melted) <- c("values", "position", "base")
levels(z_melted$base) <- c("A", "T", "C", "G")

#Plot heat map of base and corresponding LASSO coefficients
ggplot(data = z_melted, aes(x = position, y = base, fill = values)) + geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", na.value = "white", midpoint = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  geom_vline(xintercept=30.5, linetype = 2) + geom_vline(xintercept=33.5, linetype = 2) +
  ggtitle("LASSO Coefficients") + xlab("Position") + ylab("Base") + labs(fill='Values') +
  annotate("text", x = c(32), y='C', label = "PAM", angle = 90 )
  