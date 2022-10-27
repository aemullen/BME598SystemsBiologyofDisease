# Set working directory
setwd('/Users/andrewmullen/OneDrive - Arizona State University/BME 598 Systems Biology of Disease/Final Project')

library(Biobase)
library(GEOquery)
library(limma)

## Data Sets
gset1 <- getGEO("GSE52166", GSEMatrix = TRUE)
gset2 <- getGEO("GSE18323", GSEMatrix = TRUE)[[2]]
gset3 <- getGEO("GSE89292", GSEMatrix = TRUE)[[1]]
train1 = exprs(gset2)
dim(train1)

#train1 = exprs(gset)
#dim(train1)

# Parse phenotype data into classifier and potentially useful covariates
write.csv(gset2@phenoData@data,'trainingSet_gse18323_phenotypes.csv')
pheno.train1 = cbind(status = gsub('state: ','',gset2@phenoData@data$characteristics_ch1))
rownames(pheno.train1) = rownames(gset2@phenoData@data)

# Feature selection or use all?
mads = apply(train1,1,mad)
tmp = train1[order(mads,decreasing=T)[1:1000],]

# Make data object for classification
levels(as.factor(pheno.train1[,'status']))
#tmp = cbind(t(tmp),status=as.numeric(as.factor(pheno.train1[colnames(tmp),'status'])))
#tmp = cbind(t(tmp),status=as.numeric(as.factor(pheno.train1[colnames(tmp),'status']=='group: NonProtected')))
#tmp2 = as.data.frame(t(tmp))
status=as.numeric(as.factor(pheno.train1[colnames(tmp),'status']=='group: InfectivityControl'))

#Create top scoring pair classifier
library(tspair)
tsp1 = tspcalc(tmp,status)

# Plot the top-scoring pair
pdf('tb_tspair_plot.pdf')
tspplot(tsp1)
dev.off()

# Significance of top-scoring pair
out = tspsig(tmp, status, B=50,seed=42)
out$p


# Create random forest classifier
#library(randomForest)
#set.seed(42)
#tb.rf = randomForest(as.factor(status) ~ .,data=tmp,importance=T,proximity=T)
#tb.rf$importance[order(tb.rf$importance[,4],decreasing=T),]
# Get gene names for most important features
#featureData(gset)@data[rownames(tb.rf$importance[order(tb.rf$importance[,4],decreasing=T)[1:20],]),'Symbol']
#tb.rf


#write.csv(gset1@phenoData@data, 'gse52166.csv')
#write.csv(gset2@phenoData@data,'gse18323.csv')
#write.csv(gset3@phenoData@data,'gse8292.csv')
