library(rrBLUP)
library(tidyverse)
export PATH=/programs/R-4.3.3/bin:$PATH
setwd("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction")
setwd("/local/workdir/sc3339/GP")
rm(list=ls())
geno_input <- read.csv("geno_promis_bias_od_maf_filer for all groups_GWAS with maf >0.03.csv",row.names = 1)
pheno <- read.csv("Pheno(780).csv",row.names = 1)
geno_input <- geno_input[,rownames(pheno)]
identical(colnames(geno_input),rownames(pheno))

convert_genotype <- function(genotype_matrix, ploidy) {
  normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
  return(normalized_matrix)
}#Alex Sandercock code

#tranforming genotypes
geno <- convert_genotype(geno_input, 4)

set.seed(2025)
# 50 iterations (id assignment) for all grou
ids <- list()

# Loop through iterations
for (i in 1:50) {
    count_nrow <- nrow(pheno)
    # Set number of folds
    n.fold <- 5
    # Sample ids
    ids_0 <- sample(count_nrow)
    # Divide ids into folds
    ids[[i]] <- ids_0 %% n.fold
    ids[[i]][ids[[i]] == 0] <- n.fold
  }
 

save(ids,file = "ids_50random within all groups_5fold.RData")

#rrblup genomic prediction
accuracy=matrix(nrow=length(ids),ncol = 1) #prediction accuracy
pred_y <- list()
accuracy_all <- list()
file_ids <- list()
file_all <- list()

#ids_all[[2]][[1]]

for (p in 1:ncol(pheno)) {
  #entire population
  for (i in 1:length(ids)) {
    for (j in 1:n.fold) {
      yTrain <- pheno[,p]
      yTrain[ids[[i]] == j] <- NA
      gTrain <- t(geno)
      gTrain[ids[[i]] == j] <- NA
      y_model <- mixed.solve(yTrain,Z=gTrain,K=NULL,SE=FALSE,return.Hinv = FALSE)
      e=as.matrix(y_model$u)
      pred_y_GEBV=t(geno)[ids[[i]]==j,]%*%e
      pred_y[[j]] = pred_y_GEBV[,1]+as.numeric(y_model$beta)
      pred_all <- data.frame(unlist(pred_y))
      pred_all <- pred_all[rownames(pheno),]
      pred_y_GEBV_all <-  pred_all-as.numeric(y_model$beta)
    }
    #write.table(data.frame(ID=rownames(pheno_all),"predicted value"=pred_all),file = paste0("~/Desktop/test for GP/predicted value_iteration_",i,"_",names(pheno_all)[p],".csv"),row.names = F)
    #write.table(data.frame(ID=rownames(pheno_all),GEBV=pred_y_GEBV_all,"predicted value"=pred_all),file = paste0("/local/workdir/sc3339/GP/predicted values/predicted value_iteration_",i,"_",names(pheno_all)[p],".csv"),row.names = F)
    file_ids[[i]]<-  data.frame(ID=rownames(pheno),"observed value"=as.numeric(pheno[,p]),GEBV=pred_y_GEBV_all,"predicted value"=pred_all)
    accuracy[i,1]<-cor(as.numeric(pheno[,p]), pred_all, use = "complete.obs")
  }
  file_all[[p]] <- file_ids
  names(file_all)[p] <- names(pheno[p])
  accuracy_all[[p]] <- accuracy
  names(accuracy_all)[p] <- names(pheno[p])
}
save(file_all,file = "predicted values for entirePOP_5fold.RData")
save(accuracy_all,file = "accuracy for entirePOP_5fold.RData")

