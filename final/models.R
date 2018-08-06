library(e1071)
library(rpart)
library(randomForest)
#install.packages('caret', dependencies = TRUE)
#install.packages("fields")
library(caret)
library(fields)

# load the CSV file from the local directory
filename <- "/Users/JYang/Desktop/Stanford/lasso_data.csv"
dataset_orig <- read.csv(filename, header=TRUE)

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

pred_var = test_df$pos_vs_neg*test_df$avg_log2fc

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


###########Support Vector Model
svr.model <- svm(y=pred_var, x=x)
fit_test <- predict(svr.model, newX)
error <- validation$pred_val - fit_test
mse <- mean(error^2)

#install.packages("ggpubr")
library(ggpubr)

pearson_corr = cor(validation$pred_val, fit_test, method = "pearson") #"spearman"))

library(ggplot2)
df <- data.frame(cbind(validation$pred_val, fit_test))
colnames(df)<-c('Actual', 'Predicted')
df$pc <- predict(prcomp(~validation$pred_val+fit_test, df))[,1]
p <- ggplot(df, aes(Actual, Predicted, color = pc)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.1) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") + 
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("Support Vector Regression - Log Fold Change") + xlab("Actual") + ylab("Predicted") +
  annotate("text", x = c(-0.3), y=c(0.1), label = paste0("Pearson Correlation: ", "\n", pearson_corr))
p


###########Random Forest Model

rf.model <- randomForest(y=pred_var,x= x)

fit_test<-predict(rf.model, newX)
plot(validation$pred_val, fit_test, col = "blue", pch=4)
error <- validation$pred_val - fit_test
mse <- mean(error^2)

#install.packages("ggpubr")
#library(ggpubr)

pearson_corr = cor(validation$pred_val, fit_test, method = "pearson") #"spearman"))

library(ggplot2)
df <- data.frame(cbind(validation$pred_val, fit_test))
colnames(df)<-c('Actual', 'Predicted')
df$pc <- predict(prcomp(~validation$pred_val+fit_test, df))[,1]
p <- ggplot(df, aes(Actual, Predicted, color = pc)) + theme_minimal() +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = 0.1) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") + 
  theme(plot.title = element_text(hjust = 0.5), title =element_text(size=12, face='bold')) + 
  ggtitle("Random Forest Regression - Log Fold Change") + xlab("Actual") + ylab("Predicted") +
  annotate("text", x = c(-0.3), y=c(0.1), label = paste0("Pearson Correlation: ", "\n", pearson_corr))
p
