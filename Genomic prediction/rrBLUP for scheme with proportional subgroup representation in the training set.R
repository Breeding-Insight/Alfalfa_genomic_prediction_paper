library(rrBLUP)
library(tidyverse)
export PATH=/programs/R-4.3.3/bin:$PATH
setwd("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction")
setwd("/local/workdir/sc3339/GP")


geno_input <- read.csv("geno_promis_bias_od_maf_filer for all groups_GWAS with maf >0.03.csv",row.names = 1)
pheno <- read.csv("Pheno(780).csv",row.names = 1)
geno_input <- geno_input[,rownames(pheno)]  #identical order of IDs between geno and pheno
identical(colnames(geno_input),rownames(pheno))

convert_genotype <- function(genotype_matrix, ploidy) {
  normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
  return(normalized_matrix)
}#Alex Sandercock code

#tranforming genotypes
geno <- convert_genotype(geno_input, 4)


Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
#combine five groups for genotype and phenotype
pheno_group <- list()
geno_group <- list()
for (i in 1:length(Group)) {
  pheno_group[[i]] <- pheno[grep(pattern = Group[i], rownames(pheno)), ]
  geno_group[[i]] <- geno[,grep(pattern = Group[i], colnames(geno))]
  pheno_all <- do.call(rbind,pheno_group)
  geno_all <- t(do.call(cbind,geno_group))
}

#genomic prediction
set.seed(2025)
# Initialize lists to store results
cycles <- 50
file_percen <- list()
accuracy <- matrix(NA, nrow=9, ncol=length(pheno_all))
rownames(accuracy) <- paste0("train_", seq(90, 10, by = -10))
colnames(accuracy) <- names(pheno_all)

file_pheno <- list()
pheno_train_pheno <- list()

file_cycle <- list()
pheno_train_all <- list()
accuracy_cycle <- list()

pheno_test <- list()
train_90_group <- list()
train_80_group <- list()
train_70_group <- list()
train_60_group <- list()
train_50_group <- list()
train_40_group <- list()
train_30_group <- list()
train_20_group <- list()
train_10_group <- list()

for (r in 1:cycles) {
  #scheme 6
  for (i in 1:length(Group)) {
  test <- sample(1:nrow(pheno_group[[i]]), round(0.1 * nrow(pheno_group[[i]])))
  pheno_test[[i]] <- pheno_group[[i]][test, ]
  pheno_test_all <- do.call(rbind,pheno_test)
  
  train_90 <- setdiff(1:nrow(pheno_group[[i]]), test)
  train_90_group[[i]] <- pheno_group[[i]][train_90, ]
  train_90_all <- do.call(rbind,train_90_group)
  
  train_80 <- sample(train_90, round(0.8 * nrow(pheno_group[[i]])))
  train_80_group[[i]] <- pheno_group[[i]][train_80, ]
  train_80_all <- do.call(rbind,train_80_group)
    
  train_70 <- sample(train_80, round(0.7 * nrow(pheno_group[[i]])))
  train_70_group[[i]] <- pheno_group[[i]][train_70, ]
  train_70_all <- do.call(rbind,train_70_group)
  
  train_60 <- sample(train_70, round(0.6 * nrow(pheno_group[[i]])))
  train_60_group[[i]] <- pheno_group[[i]][train_60, ]
  train_60_all <- do.call(rbind,train_60_group)
  
  train_50 <- sample(train_60, round(0.5 * nrow(pheno_group[[i]])))
  train_50_group[[i]] <- pheno_group[[i]][train_50, ]
  train_50_all <- do.call(rbind,train_50_group)
  
  train_40 <- sample(train_50, round(0.4 * nrow(pheno_group[[i]])))
  train_40_group[[i]] <- pheno_group[[i]][train_40, ]
  train_40_all <- do.call(rbind,train_40_group)
  
  train_30 <- sample(train_40, round(0.3 * nrow(pheno_group[[i]])))
  train_30_group[[i]] <- pheno_group[[i]][train_30, ]
  train_30_all <- do.call(rbind,train_30_group)
  
  train_20 <- sample(train_30, round(0.2 * nrow(pheno_group[[i]])))
  train_20_group[[i]] <- pheno_group[[i]][train_20, ]
  train_20_all <- do.call(rbind,train_20_group)
  
  train_10 <- sample(train_20, round(0.1 * nrow(pheno_group[[i]])))
  train_10_group[[i]] <- pheno_group[[i]][train_10, ]
  train_10_all <- do.call(rbind,train_10_group)
}
  
percentage <- list(train_90_all,train_80_all,train_70_all,train_60_all,train_50_all,train_40_all,train_30_all,train_20_all,train_10_all)
names(percentage) <- paste0("train_", seq(90, 10, by = -10))

for (p in 1:length(pheno_all)) {
for (j in 1:length(percentage)) {
  y <- percentage[[j]][,p]
  geno_train_all <- geno_all[which(rownames(geno_all)%in%rownames(percentage[[j]])),]
  geno_train_all <-  geno_train_all[rownames(percentage[[j]]),]
  print(identical(rownames(geno_train_all),rownames(percentage[[j]])))
  y_model <- mixed.solve(y, Z=geno_train_all, K=NULL, SE=FALSE, return.Hinv=FALSE)
  e <- as.matrix(y_model$u)
  geno_test <- geno_all[which(rownames(geno_all)%in%rownames(pheno_test_all)),]
  pred_y_GEBV <- geno_test%*% e
  pred_y <- pred_y_GEBV[, 1] + as.numeric(y_model$beta)
  pred_y <- data.frame(pred_y)
  pred_y <- pred_y[rownames(pheno_test_all), ]
  pred_y_GEBV_all <- pred_y - as.numeric(y_model$beta)
  y_valid <- pheno_test_all[, p]
  
  file_percen[[j]] <- data.frame(ID=rownames(pheno_test_all), "observed value"=y_valid,GEBV=pred_y_GEBV_all,
                                 "predicted value"=pred_y)
  names(file_percen)[j] <- names(percentage)[j]
  accuracy[j, p] <- cor(pred_y, y_valid, use="complete")
  
}
 file_pheno[[p]] <- file_percen
 names(file_pheno)[p] <- names(pheno_all)[p]
 pheno_train_pheno[[p]] <- percentage
 names(pheno_train_pheno)[p] <- names(pheno_all)[p]
  
 }
  
  file_cycle[[r]] <- file_pheno
  pheno_train_all[[r]] <- pheno_train_pheno
  accuracy_cycle[[r]] <- accuracy
}


save(file_cycle,file = "file_cycle for scheme6.RData")
save(pheno_train_all,file = "pheno_train_all for scheme6.RData")
save(accuracy_cycle,file = "accuracy_cycle for scheme6.RData")

load("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme6/accuracy_cycle for scheme6.RData")
combined_data <- data.frame(do.call(rbind, accuracy_cycle))
combined_data$scheme <- rep(paste0(seq(90,10,by=-10),"%"),times=50)
# combined_data$scheme <- rep(paste0(seq(80,10,by=-10)),times=50)
combined_data_melt <- reshape::melt(combined_data,id.vars="scheme")
#write.table(combined_data,file = "combined_data for scheme4_OTTM.txt",row.names = F)
#combined_data_order <- combined_data[order(combined_data$scheme,decreasing = T),]


give.n<- function(y) {
  return( 
    data.frame(
      y =1.1*max(y),  #may need to modify this depending on your data
      label = paste(#'n = ', length(y), '\n',
        #'mean =', 
        round(mean(y), 3), '\n'
      )
    )
  )
}


png("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme6/Scheme6_3.png",height = 1500,width = 4500,res = 300)
ggplot(combined_data_melt,aes(y=combined_data_melt$value,x=combined_data_melt$scheme,fill=scheme))+
  geom_boxplot(width = 0.4)+
  facet_wrap(~factor(variable,levels = c("GH","HGT","VIG")),dir = "h",ncol = 4)+
  theme_bw()+xlab("Percentage of training set (%)")+ylab("Predictive ability")+
  #ggtitle("EURO")+
  scale_y_continuous(limits = c(0,1))+
  theme(text = element_text(size = 15),legend.position = "none")+
# stat_summary(fun.data = give.n, geom = "text", 
#              aes(label = ..label.., color = I(..color..), vjust = ..vjust..), 
#              fun = mean)+
  stat_summary(fun.data = give.n, geom = "text")+
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") 

dev.off()


