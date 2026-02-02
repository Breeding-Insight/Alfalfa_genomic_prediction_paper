library(rrBLUP)
library(tidyverse)
export PATH=/programs/R-4.3.3/bin:$PATH
setwd("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction")
setwd("/local/workdir/sc3339/GP")

pheno <- read.csv("Pheno(780).csv",row.names = 1)
geno_input <- read.csv("geno_promis_bias_od_maf_filer for all groups_GWAS with maf >0.03.csv",row.names = 1)
geno_input <- t(geno_input)
geno_input <- geno_input[rownames(pheno),] #identical order of IDs between geno and pheno
identical(rownames(geno_input),rownames(pheno))

convert_genotype <- function(genotype_matrix, ploidy) {
  normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
  return(normalized_matrix)
}#Alex Sandercock code

#tranforming genotypes
geno <- convert_genotype(geno_input, 4)

set.seed(2025)
#rrblup genomic prediction
Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
pred_y <- list()
accuracy_all <- list()
file_all <- list()
n.fold=5

for (p in 1:ncol(pheno)) {
  for (i in 1:n.fold) {
    yTrain <-  pheno[-grep(pattern = Group[i], rownames(pheno)), ][,p]
    gTrain <- geno[-grep(pattern = Group[i], rownames(geno)),]
    print(identical(rownames(pheno[-grep(pattern = Group[i], rownames(pheno)),]),rownames(gTrain)))
    y_model <- mixed.solve(yTrain,Z=gTrain,K=NULL,SE=FALSE,return.Hinv = FALSE)
    e=as.matrix(y_model$u)
    pred_y_GEBV=geno[grep(pattern = Group[i], rownames(geno)),]%*%e
    pred_y[[i]] = pred_y_GEBV[,1]+as.numeric(y_model$beta)
    pred_all <- data.frame(unlist(pred_y))
    pred_all <- pred_all[rownames(pheno),]
    pred_y_GEBV_all <-  pred_all-as.numeric(y_model$beta)
  }
  #write.table(data.frame(ID=rownames(pheno),"predicted values"=pred_all,file = paste0("~/Desktop/test for GP/predicted value_",names(pheno_all)[p],".txt"))
  file_all[[p]]<-  data.frame(ID=rownames(pheno),"observed value"=as.numeric(pheno[,p]),GEBV=pred_y_GEBV_all,"predicted value"=pred_all)
  names(file_all)[p] <- names(pheno)[p]
  accuracy_all[[p]]<-cor(as.numeric(pheno[,p]), pred_all, use = "complete.obs")
  names(accuracy_all)[p] <- names(pheno)[p]
}
save(file_all,file = "predicted values for scheme 1.rda")
save(accuracy_all,file = "accuracy for scheme1.rda")

load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme1/predicted values for scheme1.rda")
#accuracy for each group
accuracy <- list()
for (i in 1:length(Group)) {
  accuracy[[i]] <- list()
  for (j in 1:ncol(pheno)) {
    group_value <- file_all[[j]][grep(pattern = Group[i],file_all[[j]]$ID),]
    accuracy[[i]][[j]] <- cor(group_value$observed.value,group_value$predicted.value, use = "complete.obs")
    names(accuracy[[i]])[[j]]<- names(pheno)[j]
  }
  names(accuracy)[[i]]<- Group[i]
}


