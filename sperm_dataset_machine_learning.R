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

#install.packages("Boruta")
library(Boruta) # random forest feature selection
library(sigFeature) # for svm feature selection 

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
cmRf <- confusionMatrix(predRf, as.factor(dmr20$Tertile), positive = "Third")
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
cmSvm <- confusionMatrix(predSvm, as.factor(dmr20$Tertile), positive = "Third")
# 2/20 misclassified
# Accuracy : 0.9
# Kappa : 0.8
# No Information Rate : 0.5 
# P-Value [Acc > NIR] : 0.0002012 


cmTable <- function(cm, modelType) {
  if(modelType == "rf") {
    title <- c("RF confusion matrix" = 4)
  }
  if(modelType == "svm") {
    title <- c("SVM confusion matrix" = 4)
  }
  df <- data.frame(. = c("Predicted"),
                   .. = c("Third", "First"),
                   Third = c(cm$table["Third", "Third"], cm$table["First", "Third"]),
                   First = c(cm$table["Third", "First"], cm$table["First", "First"]))
  df %>%
    kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
    kable_styling(font_size = 14, full_width = F) %>%
    add_header_above(c(" ", " ", "Actual" = 2)) %>%
    add_header_above(header = title, align = "c") %>%
    column_spec(column = 1:2, bold = TRUE, color = "black") %>%
    collapse_rows(columns = 1)
}

cmKableSvm <- cmTable(cmSvm, "svm") 
cmKableSvm
cmKableRf <- cmTable(cmRf, "rf")

cmTableLegend <- function() {
  df <- data.frame(. = c("Predicted"),
                   .. = c("Positive_class", "Negative_class"),
                   Positive_class = c("True positive", "False negative"),
                   Negative_class = c("False positive", "True negative"))
  df %>%
    kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
    kable_styling(font_size = 14, full_width = F) %>%
    add_header_above(c(" ", " ", "Actual" = 2)) %>%
    add_header_above(header = c("Confusion matrix legend" = 4), align = "c") %>%
    column_spec(column = 1:2, bold = TRUE, color = "black") %>%
    collapse_rows(columns = 1)
}

cmLegend <- cmTableLegend()

#save_kable(cmKableSvm, "cmKableSvm.pdf")
#save_kable(cmKableRf, "cmKableRf.pdf")

resRf <- c(format(cmRf$overall["Accuracy"], 2),
           0.5, 
           0.0012880, 
           cmRf$overall["Kappa"])
resSvm <- c(format(cmSvm$overall["Accuracy"], 2),
            0.5, 
            0.0002012, 
            cmSvm$overall["Kappa"])

res <- data.frame(RF = resRf, SVM = resSvm)
rownames(res) <- c("Accuracy", "No Information Rate (NIR)", "P-value [Acc > NIR]", "Kappa")
#res
resKable <- res %>%
  kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
  kable_styling(font_size = 14, full_width = F) %>%
  add_header_above(header = c("Summary of model results" = 3), align = "c") %>%
  column_spec(column = 1, bold = TRUE)


# RF variable importance
varImpRf <- varImp(object = modelRf, scale = FALSE)
varImpRfList <- varImpRf$importance
rownames(varImpRfList) <- rownames(varImpRfList) %>% str_remove_all("`")
varImpRfList <- data.frame(DMR = rownames(varImpRfList),
                           Variable_Importance_Measure = varImpRfList$Overall) %>%
  arrange(desc(Variable_Importance_Measure)) %>%
  add_column(Ranking = 1:nrow(varImpRfList), .before = 1)
#varImpRfList 
  
# SVM variable importance 
# filterVarImp returns area under ROC curve
varImpSvm <- filterVarImp(x = dmr32[,-1], y = dmr32$Tertile) 
varImpSvmList <- data.frame(DMR = rownames(varImpSvm),
                            Variable_Importance_Measure = varImpSvm$First) %>%
  arrange(desc(Variable_Importance_Measure)) %>%
  add_column(Ranking = 1:nrow(varImpSvm), .before = 1)
#varImpSvmList


varImpTable <- function(varImpList, modelType, colNum) {
  if(modelType == "rf") {
    title <- c("RF (varImp) variable importance list" = colNum)
  }
  if(modelType == "svm") {
    title <- c("SVM (filterVarImp) variable importance list" = colNum)
  }
  if(modelType == "boruta") {
    title <- c("RF (Boruta) variable importance list" = colNum)
  }
  if(modelType == "sigFeature") {
    title <- c("SVM (sigFeature) variable importance list" = colNum)
  }
  varImpList %>%
    kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
    kable_styling(font_size = 14, full_width = F) %>%
    add_header_above(header = title, align = "c") %>%
    column_spec(column = 1:colNum, color = "black") 
}

rfKable <- varImpTable(varImpRfList, "rf", 3)
svmKable <- varImpTable(varImpSvmList, "svm", 3)

#save_kable(varImpRfList, "varImpRfList.pdf")
#save_kable(varImpSvmList, "varImpSvmList.pdf")


# https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=2ahUKEwj4i4Gy4rTiAhVM6Z8KHZiKAnMQFjAAegQIBhAC&url=https%3A%2F%2Fwww.jstatsoft.org%2Farticle%2Fview%2Fv036i11%2Fv36i11.pdf&usg=AOvVaw3tyiHN0BCe2fkkAA6xEVDE
# https://academic.oup.com/bib/article/20/2/492/4554516
# http://r-statistics.co/Variable-Selection-and-Importance-With-R.html
# Boruta RF Variable Importance

set.seed(seed)
boruta.train <- Boruta(Tertile ~ ., data = dmr32, doTrace=2)  
#boruta.train
#plot(boruta.train)
boruta.train.stats <- attStats(boruta.train)

final.boruta <- TentativeRoughFix(boruta.train)
#final.boruta

borutaVars <- getSelectedAttributes(final.boruta)
#borutaVars

# borutaVars but ranked - only final 17
boruta.df <- attStats(final.boruta)
borutaRF <- boruta.df[which(boruta.df$decision == "Confirmed"), ]
rownames(borutaRF) <- rownames(borutaRF) %>% str_remove_all("`")
borutaVarsRf <- data.frame(DMR = rownames(borutaRF),
                           meanImp = borutaRF$meanImp) %>% 
  arrange(desc(meanImp)) %>%
  add_column(Ranking = 1:nrow(borutaRF), .before = 1)
#borutaVarsRf

# ranked boruta for all
boruta.train.stats <- attStats(boruta.train)
rownames(boruta.train.stats) <- rownames(boruta.train.stats) %>% str_remove_all("`")
borutaAllVarsRf <- data.frame(DMR = rownames(boruta.train.stats),
                           meanImp = boruta.train.stats$meanImp) %>% 
  arrange(desc(meanImp)) %>%
  add_column(Ranking = 1:nrow(boruta.train.stats), .before = 1)
#borutaAllVarsRf
borutaKable <- varImpTable(borutaAllVarsRf, "boruta", 3)
save(borutaAllVarsRf, file = "borutaAllVarsRf.RData")

# includes "Tentative" + everything in borutaRF
# boruta.sig <- names(boruta.train$finalDecision[boruta.train$finalDecision %in% c("Confirmed", "Tentative")])  # collect Confirmed and Tentative variables
# boruta.sig

# sigFeature SVM Variable Importance
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# 
# BiocManager::install("sigFeature")
# https://www.bioconductor.org/packages/devel/bioc/vignettes/sigFeature/inst/doc/vignettes.pdf
# https://rdrr.io/bioc/sigFeature/man/sigFeature.html

dmr32Matrix <- as.matrix(dmr32[,-1])
sigFeatureRanking <- sigFeature(dmr32Matrix, dmr32$Tertile)
sigFeatureDmrs <- colnames(dmr32Matrix[, sigFeatureRanking]) %>% str_remove_all("`")
sigfeatureVarsSvm <- data.frame(Ranking = 1:length(sigFeatureDmrs),
                                DMR = sigFeatureDmrs)
#sigfeatureVarsSvm
sigfeatureKable <- varImpTable(sigfeatureVarsSvm, "sigFeature", 2)
save(sigfeatureVarsSvm, file = "sigfeatureVarsSvm.RData")

#sigfeatureKable
varsDf <- data.frame(Ranking = 1:261,
                     DMR_varImp = varImpRfList$DMR,
                     DMR_filterVarImp = varImpSvmList$DMR,
                     DMR_Boruta = borutaAllVarsRf$DMR,
                     DMR_sigFeature = sigfeatureVarsSvm$DMR)
#varsDf
varDf20 <- varsDf[1:20,]
#varDf20

# top20Overlap_var <- intersect(varDf20$DMR_varImp, varDf20$DMR_filterVarImp)
# top20Overlap_var # 10/20 overlap

# 15/20 DMRs overlap in top 20 DMR list for Boruta and sigFeature
top20Overlap <- data.frame(Number = 1:15,
                           DMR = intersect(varDf20$DMR_Boruta, varDf20$DMR_sigFeature))
#top20Overlap 

# Data frame for top 20 DMRs, overlapping ones have * character appended to DMR name
top20 <- data.frame(Ranking = 1:20,
                    Boruta = as.character(varDf20$DMR_Boruta), 
                    sigFeature = as.character(varDf20$DMR_sigFeature)) %>%
  mutate_if(is.factor, as.character) 
top20 <- top20 %>%
  mutate(Boruta = replace(Boruta, 
                          which(top20$Boruta %in% top20Overlap$DMR == TRUE),
                          paste(Boruta[which(top20$Boruta %in% top20Overlap$DMR == TRUE)], 
                                "*",
                                sep=""))) %>%
  mutate(sigFeature = replace(sigFeature, 
                          which(top20$sigFeature %in% top20Overlap$DMR == TRUE),
                          paste(sigFeature[which(top20$sigFeature %in% top20Overlap$DMR == TRUE)], 
                                "*",
                                sep="")))

# Kable table for top 20 DMRs
top20Kable <- top20 %>%
  kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
  kable_styling(font_size = 14, full_width = F) %>%
  add_header_above(header = c("Top 20 DMR predictors" = 3), align = "c") %>%
  column_spec(column = 1:3, color = "black") 
#top20Kable

# Kable table for 15 overlapping DMRs in top 20 DMRs
top20OverlapKable <- top20Overlap %>%
  kable(table.attr = "style = \"color: black; font-family: Calibri, sans-serif\"") %>%
  kable_styling(font_size = 14, full_width = F) %>%
  add_header_above(header = c("Overlapping top 20 DMR predictors" = 2), align = "c") %>%
  column_spec(column = 1:2, color = "black") 
#top20OverlapKable


cmRf
cmSvm

cmLegend
resKable

cmKableSvm
cmKableRf

borutaKable
sigfeatureKable

top20Kable
top20OverlapKable

save_kable(cmLegend, "cmLegend.pdf")
save_kable(cmKableSvm, "cmKableSvm.pdf")
save_kable(cmKableRf, "cmKableRf.pdf")

load("top20Kable.RData")
save_kable(top20Kable, "top20Kable.pdf")

load("top20OverlapKable.RData")
save_kable(top20OverlapKable, "top20OverlapKable.pdf")

load("borutaKable.RData")
save_kable(borutaKable, "borutaKable.pdf")

load("sigfeatureKable.RData")
save_kable(sigfeatureKable, "sigfeatureKable.pdf")

write.csv(sigfeatureVarsSvm, "sigfeatureVarsSvm.csv")
write.csv(borutaAllVarsRf, "borutaAllVarsRf.csv")

# Data from dmr32 for sigfeature list, boruta list, top 15 list
svmDmrs <- sigfeatureVarsSvm$DMR
rfDmrs <- borutaAllVarsRf$DMR
overlap15Dmrs <- top20Overlap$DMR

save(svmDmrs, file = "svmDmrs.RData")
save(rfDmrs, file = "rfDmrs.RData")
save(top20Overlap, file = "top20Overlap.RData")

svmDmrs %in% dmr32

dmr32NoStr <- colnames(dmr32)  %>% str_remove_all("`")


svmDmrs %in% dmr32


















