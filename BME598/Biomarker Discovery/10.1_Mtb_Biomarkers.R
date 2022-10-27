##########################################################
## Systems Biology of Disease: MTB Biomarker Discovery  ##
## Tutorial from:  Satija Lab                           ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Adapted by: Plaisier Lab                            ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##########################################################

# Set working directory
setwd('/Users/andrewmullen/OneDrive - Arizona State University/BME 598 Systems Biology of Disease/Biomarker Discovery')

#################################
### Load Latent vs. Active TB ###
#################################
library(Biobase)
library(GEOquery)
library(limma)

#####################################################
### The Blood Transcriptome of Pulmonary Diseases ###
### (PMID = 23940611)                             ###
#####################################################
## Training set
gset <- getGEO("GSE42830", GSEMatrix =TRUE)[[1]]
#gset <- getGEO('GSE33940')
train1 = exprs(gset)
dim(train1)

# Parse phenotype data into classifier and potentially useful covariates
write.csv(gset@phenoData@data,'trainingSet_gse42830_phenotypes.csv')
diseaseState = c()
for(i in gset@phenoData@data[,'title']) {
  diseaseState = c(diseaseState, strsplit(i,' ')[[1]][1])
}
pheno.train1 = cbind(status = diseaseState,
                     gender = gsub('gender: ','',gset@phenoData@data$characteristics_ch1),
                     ethnicity = gsub('ethnicity: ','',gset@phenoData@data$characteristics_ch1.1),
                     disease = gsub('disease state: ','',gset@phenoData@data$characteristics_ch1.2))
rownames(pheno.train1) = rownames(gset@phenoData@data)

# Feature selection or use all?
mads = apply(train1,1,mad)
tmp = train1[order(mads,decreasing=T)[1:1000],]

# Make data object for classification
levels(as.factor(pheno.train1[,'status']))
#tmp = cbind(t(tmp),status=as.numeric(as.factor(pheno.train1[colnames(tmp),'status'])))
tmp = cbind(t(tmp),status=as.numeric(as.factor(pheno.train1[colnames(tmp),'status']=='TB')))

# Create random forest classifier
library(randomForest)
set.seed(42)
tb.rf = randomForest(as.factor(status) ~ .,data=tmp,importance=T,proximity=T)
tb.rf$importance[order(tb.rf$importance[,4],decreasing=T),]
# Get gene names for most important features
featureData(gset)@data[rownames(tb.rf$importance[order(tb.rf$importance[,4],decreasing=T)[1:20],]),'Symbol']
tb.rf


# Plot variable importance
pdf('tb_variable_importance_gini.pdf')
varImpPlot(tb.rf, type=2, main='Vairable Importance for Top 10 Predictors\n(Mean decrease in Gini node impurity)')
dev.off()

# Calculate ROC AUC
library(ROCR)
predictions = as.vector(tb.rf$votes[,2])
pred1 = prediction(predictions, tmp[,'status'])
auc1 = performance(pred1, 'auc')
auc1@y.values[[1]]

# Calculate predictive value of classifier
confusion1 = tb.rf$confusion[c(2,1),c(2,1)]
ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
ppv1
npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
npv1

# Plot ROC AUC for classifier
pdf('tb_ROC_curve.pdf')
perf1 = performance(pred1, 'tpr','fpr')
plot(perf1, main='ROC Curve miRNA Predictors', col='red', lwd=2)
text(0.5,0.5,paste('AUC = ',format(auc1@y.values[[1]],digits=5,scientific=FALSE),'\nPPV = ',format(ppv1,digits=5,scientific=FALSE),'\nNPV = ',format(npv1,digits=5,scientific=FALSE)),cex=2)
dev.off()


## Testing set
#gset.test1 <- getGEO("GSE42826", GSEMatrix =TRUE)[[1]]
gset.test1 <- getGEO(filename='data/GSE42826_series_matrix.txt.gz')
test1 = exprs(gset.test1)
dim(test1)
write.csv(gset.test1@phenoData@data,'testingSet_gse42826_phenotypes.csv')
diseaseState = c()
for(i in gset.test1@phenoData@data[,'title']) {
  diseaseState = c(diseaseState, strsplit(i,' ')[[1]][1])
}
pheno.test1 = cbind(status = diseaseState,
                    gender = gsub('gender: ','',gset.test1@phenoData@data$characteristics_ch1),
                    ethnicity = gsub('ethnicity: ','',gset.test1@phenoData@data$characteristics_ch1.1),
                    disease = gsub('disease state: ','',gset.test1@phenoData@data$characteristics_ch1.2))
rownames(pheno.test1) = rownames(gset.test1@phenoData@data)

# Feature selection
tmp.test1 = test1[order(mads,decreasing=T)[1:1000],]
test1.tb = as.numeric(as.factor(pheno.test1[colnames(tmp.test1),'status']=='TB'))

# Make data object for classification
tmp.test1 = cbind(t(tmp.test1),status=test1.tb)

# Getting predictions for testing hold-out data
test1.resp = predict(tb.rf, tmp.test1, type='response')
test1.vote = predict(tb.rf, tmp.test1, type='vote')

# Create confusion matrix
confusion1 = table(data.frame(cbind(Actual=test1.tb,Predicted=test1.resp)))
rownames(confusion1) = c('Not-TB','TB')
colnames(confusion1) = c('Not-TB','TB')
ro1 = c('TB','Not-TB')
confusion1 = confusion1[ro1,ro1]
confusion1
err1 = ((confusion1[1,2]+confusion1[2,1])/sum(confusion1))
err1

# Calculate ROC AUC
test1.predictions = as.vector(test1.vote[,2])
test1.pred1 = prediction(test1.predictions, test1.tb)
test1.auc1 = performance(test1.pred1, 'auc')
test1.auc1@y.values[[1]]

# Calculate predictive value of classifier
ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
ppv1
npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
npv1

# Plot ROC AUC for classifier
pdf('testing_TB_ROC_curve.pdf')
test1.perf1 = performance(test1.pred1, 'tpr','fpr')
plot(test1.perf1, main='ROC Curve', col='blue', lwd=2)
text(0.5,0.5,paste('AUC = ',format(test1.auc1@y.values[[1]],digits=5,scientific=FALSE),'\nPPV = ',format(ppv1,digits=5,scientific=FALSE),'\nNPV = ',format(npv1,digits=5,scientific=FALSE)),cex=2)
dev.off()


## Validation set
#gset.val1 <- getGEO("GSE42825", GSEMatrix =TRUE)[[1]]
gset.val1 <- getGEO(filename='data/GSE42825_series_matrix.txt.gz')
val1 = exprs(gset.val1)
dim(val1)
write.csv(gset.val1@phenoData@data,'validationSet_gse42825_phenotypes.csv')
diseaseState = c()
for(i in gset.val1@phenoData@data[,'title']) {
  diseaseState = c(diseaseState, strsplit(i,' ')[[1]][1])
}
pheno.val1 = cbind(status = diseaseState,
                   gender = gsub('gender: ','',gset.val1@phenoData@data$characteristics_ch1),
                   ethnicity = gsub('ethnicity: ','',gset.val1@phenoData@data$characteristics_ch1.1),
                   disease = gsub('disease state: ','',gset.val1@phenoData@data$characteristics_ch1.2))
rownames(pheno.val1) = rownames(gset.val1@phenoData@data)


# Feature selection
tmp.val1 = val1[order(mads,decreasing=T)[1:1000],]
val1.tb = as.numeric(as.factor(pheno.val1[colnames(tmp.val1),'status']=='TB'))

# Make data object for classification
tmp.val1 = cbind(t(tmp.val1),status=val1.tb)

# Getting predictions for testing hold-out data
val1.resp = predict(tb.rf, tmp.val1, type='response')
val1.vote = predict(tb.rf, tmp.val1, type='vote')

# Create confusion matrix
confusion1 = table(data.frame(cbind(Actual=val1.tb,Predicted=val1.resp)))
rownames(confusion1) = c('Not-TB','TB')
colnames(confusion1) = c('Not-TB','TB')
ro1 = c('TB','Not-TB')
confusion1 = confusion1[ro1,ro1]
confusion1
err1 = ((confusion1[1,2]+confusion1[2,1])/sum(confusion1))
err1

# Calculate ROC AUC
val1.predictions = as.vector(val1.vote[,2])
val1.pred1 = prediction(val1.predictions, val1.tb)
val1.auc1 = performance(val1.pred1, 'auc')
val1.auc1@y.values[[1]]

# Calculate predictive value of classifier
ppv1 = ((confusion1[1,1])/(confusion1[1,1]+confusion1[2,1]))
ppv1
npv1 = ((confusion1[2,2])/(confusion1[1,2]+confusion1[2,2]))
npv1

# Plot ROC AUC for classifier
pdf('validation_TB_ROC_curve.pdf')
val1.perf1 = performance(val1.pred1, 'tpr','fpr')
plot(val1.perf1, main='ROC Curve', col='blue', lwd=2)
text(0.5,0.5,paste('AUC = ',format(val1.auc1@y.values[[1]],digits=5,scientific=FALSE),'\nPPV = ',format(ppv1,digits=5,scientific=FALSE),'\nNPV = ',format(npv1,digits=5,scientific=FALSE)),cex=2)
dev.off()


# Functional Gene Expression Signatures in TB (PMID = 22046420)
#gset <- getGEO("GSE28623", GSEMatrix =TRUE)[[1]]
gset <- getGEO(filename='data/GSE28623_series_matrix.txt.gz')
fgs1 = exprs(gset)
dim(fgs1)
write.csv(gset@phenoData@data,'gse28623_phenotypes.csv')

# Common Patterns in TB (PMID = 22547807)
#gset <- getGEO("GSE34608", GSEMatrix =TRUE)[[1]]
gset <- getGEO(filename='data/GSE34608-GPL6480_series_matrix.txt.gz') # mRNA
cp1 = exprs(gset)
dim(cp1)
write.csv(gset@phenoData@data,'gse34608_phenotypes.csv')

# Early Detection of Tuberculosis Treatment Response (PMID = 23056259)
#gset <- getGEO("GSE40553", GSEMatrix =TRUE)[[1]]
gset <- getGEO(filename='data/GSE40553_series_matrix.txt.gz')
dr1 = exprs(gset)
dim(dr1)
write.csv(gset@phenoData@data,'gse40553_phenotypes.csv')

# Diagnosis of childhood tuberculosis and host RNA expression in Africa. (PMID = 24785206)
#gset <- getGEO('GSE39940',GSEMatrix=TRUE)
gset <- getGEO(filename='data/GSE39940_series_matrix.txt.gz')
chtb1 = exprs(gset)
dim(chtb1)
write.csv(gset@phenoData@data,'gse39940_phenotypes.csv')

# Validation study (PMID = 24785206)
#gset <- getGEO('GSE39939',GSEMatrix=TRUE)
gset <- getGEO(filename='data/GSE39939_series_matrix.txt.gz')
chtb1 = exprs(gset)
dim(chtb1)
write.csv(gset@phenoData@data,'gse39939_phenotypes.csv')

# Detection of tuberculosis in HIV-infected and -uninfected African adults using whole blood RNA expression signatures: a case-control study. (PMID = 24167453)
#gset <- getGEO('GSE37250',GSEMatrix=TRUE)
gset <- getGEO(filename='data/GSE37250_series_matrix.txt.gz')
adtb1 = exprs(gset)
dim(adtb1)
write.csv(gset@phenoData@data,'gse37250_phenotypes.csv')


