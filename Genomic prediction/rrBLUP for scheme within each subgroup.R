library(rrBLUP)
library(tidyverse)
export PATH=/programs/R-4.3.3/bin:$PATH
setwd("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction")
setwd("/local/workdir/sc3339/GP")

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
# 50 iterations (id assignment) for all group
Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
ids <- list()
ids_all <- list()

# Loop through iterations
for (i in 1:50) {
  # Loop through groups
  for (j in 1:5) {
    # Count rows for the group
    count_nrow <- nrow(pheno[grep(pattern = Group[j], rownames(pheno)), ])
    # Set number of folds
    n.fold <- 5
    # Sample ids
    ids_0 <- sample(count_nrow)
    # Divide ids into folds
    ids[[j]] <- ids_0 %% n.fold
    ids[[j]][ids[[j]] == 0] <- n.fold
  }
  # Store ids for the iteration
  ids_all[[i]] <- ids
}

save(ids_all,file = "ids_50random for all groups.RData")

#combine five groups for genotype and phenotype
pheno_group <- list()
geno_group <- list()
for (i in 1:length(Group)) {
  pheno_group[[i]] <- pheno[grep(pattern = Group[i], rownames(pheno)), ]
  geno_group[[i]] <- geno[,grep(pattern = Group[i], colnames(geno))]
  pheno_all <- do.call(rbind,pheno_group)
  geno_all <- t(do.call(cbind,geno_group))
}


#rrblup genomic prediction
accuracy=matrix(nrow=length(ids_all),ncol = 1) #prediction accuracy
pred_y <- list()
accuracy_all <- list()
file_ids <- list()
file_all <- list()

#ids_all[[2]][[1]]

for (p in 1:ncol(pheno_group[[1]])) {
  for (i in 1:length(ids_all)) {
    for (j in 1:n.fold) {
      yTrain <- pheno_group[[1]][,p]
      yTrain[ids_all[[i]][[1]] == j] <- NA
      gTrain <- t(geno_group[[1]])
      gTrain[ids_all[[i]][[1]] == j] <- NA
      y_model <- mixed.solve(yTrain,Z=gTrain,K=NULL,SE=FALSE,return.Hinv = FALSE)
      e=as.matrix(y_model$u)
      pred_y_GEBV=t(geno_group[[1]])[ids_all[[i]][[1]]==j,]%*%e
      pred_y[[j]] = pred_y_GEBV[,1]+as.numeric(y_model$beta)
      pred_all <- data.frame(unlist(pred_y))
      pred_all <- pred_all[rownames(pheno_group[[1]]),]
      pred_y_GEBV_all <-  pred_all-as.numeric(y_model$beta)
    }
    #write.table(data.frame(ID=rownames(pheno_all),"predicted value"=pred_all),file = paste0("~/Desktop/test for GP/predicted value_iteration_",i,"_",names(pheno_all)[p],".csv"),row.names = F)
    #write.table(data.frame(ID=rownames(pheno_all),GEBV=pred_y_GEBV_all,"predicted value"=pred_all),file = paste0("/local/workdir/sc3339/GP/predicted values/predicted value_iteration_",i,"_",names(pheno_all)[p],".csv"),row.names = F)
    file_ids[[i]]<-  data.frame(ID=rownames(pheno_group[[1]]),"observed value"=as.numeric(pheno_group[[1]][,p]),GEBV=pred_y_GEBV_all,"predicted value"=pred_all)
    accuracy[i,1]<-cor(as.numeric(pheno_group[[1]][,p]), pred_all, use = "complete.obs")
  }
  file_all[[p]] <- file_ids
  names(file_all)[p] <- names( pheno_group[[1]])[p]
  accuracy_all[[p]] <- accuracy
  names(accuracy_all)[p] <- names(pheno_group[[1]])[p]
}

Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
save(file_all,file = "predicted values for OTTM.RData")
save(accuracy_all,file = "accuracy for OTTM.RData")

load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/accuracy for check.RData")
load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/accuracy for OTTM.RData")
load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/accuracy for SIBR.RData")
load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/accuracy for CASIA.RData")
load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/accuracy for EURO.RData")
d_check <- data.frame(accuracy_all)
d_OTTM <- data.frame(accuracy_all)
d_SIBR <- data.frame(accuracy_all)
d_CASIA <- data.frame(accuracy_all)
d_EURO <- data.frame(accuracy_all)

d_check$ID <- "Check"
d_OTTM$ID <- "OTTM"
d_SIBR$ID <- "SIBR"
d_CASIA$ID <- "CASIA"
d_EURO$ID <- "EURO"
df_melt <- rbind.data.frame(d_check,d_OTTM,d_SIBR,d_CASIA,d_EURO)
combined_data_melt <- reshape::melt(df_melt,id.vars="ID")

give.n <- function(x) {
  n <- paste(round(mean(x),3))
  value <- as.numeric(paste(round(mean(x),3)))
  color <- ifelse(value < 0, "black", "black")
  vjust <- ifelse(value < 0, 1.6, -1.6)
  return(data.frame(y = value, label = n, color = color,vjust=vjust))
}

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

combined_data_melt$ID <- factor(combined_data_melt$ID,levels = c("OTTM","SIBR","CASIA","EURO","Check"))
library(ggplot2)
png("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme within each subgroup/within subgroup.png",height = 1500,width = 4500,res = 300)
ggplot(combined_data_melt,aes(combined_data_melt$value,x=combined_data_melt$ID,fill=ID))+
  geom_boxplot(width = 0.4)+
  facet_wrap(~factor(variable,levels = c("GH","HGT","VIG")))+
  scale_y_continuous(limits = c(-0.5,0.7))+
  theme_bw()+xlab("")+ylab("Predictive ability")+
  theme(text = element_text(size = 15),legend.position = "none")+
 # stat_summary(fun.data = give.n, geom = "text", 
 #            aes(label = ..label.., color = I(..color..), vjust = ..vjust..), 
 #           fun = mean)+
  stat_summary(fun.data = give.n, geom = "text")+
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") 

dev.off()


## combined ability from the entire pop

load("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme for entire pop/accuracy for entirePOP_5fold.RData")
combined_data <- data.frame(value=do.call(rbind, accuracy_all))
combined_data$Trait <- rep(c("GH","HGT","VIG"),each=50)
combined_data$ID <- "All"
names(combined_data)[2] <- "variable"
combined_data <- combined_data[,c(3,2,1)]
All_combine <- rbind.data.frame(combined_data_melt,combined_data)

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


All_combine$ID <- factor(All_combine$ID,levels = c("OTTM","SIBR","CASIA","EURO","Check","All"))
library(ggplot2)
png("/Users/sc3339-admin/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/within subgroup&entirePop_5fold.png",height = 1500,width = 4500,res = 300)
ggplot(All_combine,aes(All_combine$value,x=All_combine$ID,fill=ID))+
  geom_boxplot(width = 0.4)+
  facet_wrap(~factor(variable,levels = c("GH","HGT","VIG")))+
  scale_y_continuous(limits = c(-0.5,0.85))+
  theme_bw()+xlab("")+ylab("Predictive ability")+
  theme(text = element_text(size = 15),legend.position = "none")+
  # stat_summary(fun.data = give.n, geom = "text", 
  #            aes(label = ..label.., color = I(..color..), vjust = ..vjust..), 
  #           fun = mean)+
  stat_summary(fun.data = give.n, geom = "text")+
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") 

dev.off()

