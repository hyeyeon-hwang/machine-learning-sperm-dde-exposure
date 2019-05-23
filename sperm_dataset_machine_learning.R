# sperm ml
# dde exposure ml
# sperm dde exposure ml
# endocrine disruptor a risk for autism?
library(knitr)
library(kableExtra)

library(tidyverse)
library(caret)
library(e1071) # for train
library(randomForest) # for RF model
library(kernlab) # for SVM model

# Best machine learning results obtained by:
# 1. Using full dataset of 52 DMRs that is NOT batch corrected: 
#         'dmr52Full' = "52_sig_individual_smoothed_DMR_methylation.txt"
# 2. Splitting 'dmr52Full' into training and testing sets:
#         training set: 32 DMRs 
#         testing set: 20 DMRs
# 3. Training RF / SVM model with training set and 5-fold cross-validation
# 4. Predicting on the testing set using trained RF / SVM model 

# generated datasets using DMRichR: https://github.com/ben-laufer/DMRichR
dmr52Full <- read.delim("../machine-learning-sperm/52_sig_individual_smoothed_DMR_methylation.txt")
# generated datasets of 32 and 20 DMRs to extract sample names later
dmr32Full <- read.delim("../machine-learning-sperm/merged_32_individual_dmr.txt")
dmr20Full <- read.delim("../machine-learning-sperm/merged_20_individual_dmr.txt")

# sample info
sampleInfo52 <- read.csv("../machine-learning-sperm/sample_info.csv", fileEncoding = "UTF-8-BOM") 
# extracted sample names for 32 and 20 DMRs
sampleName32 <- colnames(dmr32Full)[-(1:3)]
sampleName20 <- colnames(dmr20Full)[-(1:3)]

# prepare / wrangle data in order to use caret() to create machine learning models
prepData <- function(dmrFull, sampleInfo) {
  # exclude columns not in first 3 (seqnames, start, end) or last 52 (samples)
  numCols <- ncol(dmrFull)
  numSamples <- 52
  excludeStart <- 4 # width column
  excludeEnd <- numCols - numSamples # percentDifference column

  # transpose data
  transposedData <- dmrFull %>%
    as.tibble() %>% 
    select(-(excludeStart:excludeEnd)) %>%
    unite(seqId1, seqnames, start, sep = ":") %>%
    unite(seqId, seqId1, end, sep = "-") %>%
    gather(SampleID, values, -seqId) %>% 
    spread(seqId, values) 
  # add 'tertile' column to transposed dataset
  addTertile <- transposedData %>%
    add_column(Tertile = as.factor(sampleInfo$Tertile[match(transposedData$SampleID, sampleInfo$Name)]), .after = 1)
  
  return(addTertile)
}

dmr52 <- prepData(dmr52Full, sampleInfo52) # dim(dmr52): 52 x 263
# training set
dmr32 <- dmr52[which(dmr52$SampleID %in% sampleName32), -1] # exclude 'SampleID' 1st column, dim(dmr32): 32 x 262
# testing set
dmr20 <- dmr52[which(dmr52$SampleID %in% sampleName20), -1] # exclude 'SampleID' 1st column, dim(dmr20): 20 x 262


# Random Forest (RF) and Support Vector Machine (SVM) models --------------

seed <- 9999

# 5-fold cross-validation
trainControl <- trainControl(method = "cv", 
                          number = 5,
                          verboseIter = TRUE,
                          returnResamp = "final",
                          savePredictions = "final", 
                          classProbs = TRUE)

# Random Forest (RF) model function
fitModelRf <- function(dmrTrainData) {
  set.seed(seed)
  modelRf <- train(Tertile ~ .,
                   data = dmrTrainData, 
                   method = "rf",
                   trControl = trainControl)
  return(modelRf)
}


# Support Vector Machine (SVM) model function
fitModelSvm <- function(dmrTrainData) {
  #dmrTrainData <- dmr32

  set.seed(seed)
  modelSvm <- train(Tertile ~ .,
                    data = dmrTrainData, # exclude 'SampleID' column
                    method = "svmLinear",
                    trControl = trainControl)
  return(modelSvm)
}



# Build RF model with training data 'dmr32'
modelRf <- fitModelRf(dmr32)
# Predict on testing data 'dmr20' with RF model
predRf <- predict(modelRf, dmr20)
# Confusion matrix of prediction results
cmRf <- confusionMatrix(predRf, as.factor(dmr20$Tertile))
# 3/20 misclassified
# Accuracy : 0.85 
# Kappa : 0.7
# No Information Rate : 0.5 
# P-Value [Acc > NIR] : 0.001288


# Build SVM model with training data 'dmr32'
modelSvm <- fitModelSvm(dmr32)
# Predict on testing data 'dmr20' with SVM model
predSvm <- predict(modelSvm, dmr20)
# Confusion matrix of prediction results
cmSvm <- confusionMatrix(predSvm, as.factor(dmr20$Tertile))
# 2/20 misclassified
# Accuracy : 0.9
# Kappa : 0.8
# No Information Rate : 0.5 
# P-Value [Acc > NIR] : 0.0002012 


cmTable <- function(cm, modelType) {
  if(modelType == "rf") {
    title <- c("RF 5-fold cross-validated confusion matrix" = 4)
  }
  if(modelType == "svm") {
    title <- c("SVM 5-fold cross-validated confusion matrix" = 4)
  }
  df <- data.frame(. = c("Predicted"),
                   .. = c("First", "Third"),
                   First = c(cm$table["First", "First"], cm$table["Third", "First"]),
                   Third = c(cm$table["First", "Third"], cm$table["Third", "Third"]))
  df %>%
    kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
    kable_styling(font_size = 14, full_width = F) %>%
    add_header_above(c(" ", " ", "Actual" = 2)) %>%
    add_header_above(header = title, align = "c") %>%
    column_spec(column = 1:2, bold = TRUE, color = "black") %>%
    collapse_rows(columns = 1)
}

cmKableSvm <- cmTable(cmSvm, "svm") 
cmKableRf <- cmTable(cmRf, "rf")

save_kable(cmKableSvm, "cmKableSvm.pdf")
save_kable(cmKableRf, "cmKableRf.pdf")

# RF variable importance
varImpRf <- varImp(object = modelRf)
varImpRfList <- varImpRf$importance
rownames(varImpRfList) <- rownames(varImpRfList) %>% str_remove_all("`")
varImpRfList <- data.frame(DMR = rownames(varImpRfList),
                           Variable_Importance_Measure = varImpRfList$Overall) %>%
  arrange(desc(Variable_Importance_Measure)) %>%
  add_column(Ranking = 1:nrow(varImpRfList), .before = 1)
varImpRfList 
  
# SVM variable importance 
# filterVarImp returns area under ROC curve
varImpSvm <- filterVarImp(x = dmr32[,-1], y = dmr32$Tertile) 
varImpSvmList <- data.frame(DMR = rownames(varImpSvm),
                            Variable_Importance_Measure = varImpSvm$First) %>%
  arrange(desc(Variable_Importance_Measure)) %>%
  add_column(Ranking = 1:nrow(varImpSvmList), .before = 1)
varImpSvmList


varImpTable <- function(varImpList, modelType) {
  if(modelType == "rf") {
    title <- c("RF Variable Importance List" = 3)
  }
  if(modelType == "svm") {
    title <- c("SVM Variable Importance List" = 3)
  }
  varImpList %>%
    kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
    kable_styling(font_size = 14, full_width = F) %>%
    add_header_above(header = title, align = "c") %>%
    column_spec(column = 1:3, color = "black") 
}

varImpTable(varImpRfList, "rf")
varImpTable(varImpSvmList, "svm")

save_kable(varImpRfList, "varImpRfList.pdf")
save_kable(varImpSvmList, "varImpSvmList.pdf")
