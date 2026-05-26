library(rrBLUP)
library(tidyverse)
export PATH=/programs/R-4.3.3/bin:$PATH
setwd("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction")
setwd("/local/workdir/sc3339/GP")

## Genomic prediction: 2016 to predict 2017 and 2018
geno_input <- read.csv("geno_promis_bias_od_maf_filer for all groups_GWAS with maf >0.03.csv",row.names = 1)
pheno <- read.csv("VIG_2016_2017_ar1_BLUE(780).csv",row.names = 1)
pheno_16 <- read.csv("VIG_2016_ar1_BLUE(780).csv",row.names = 1)
identical(rownames(pheno_16),rownames(pheno))
geno_input <- geno_input[,rownames(pheno)]  #identical order of IDs between geno and pheno
identical(colnames(geno_input),rownames(pheno))

convert_genotype <- function(genotype_matrix, ploidy) {
  normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
  return(normalized_matrix)
}#Alex Sandercock code

#tranforming genotypes
geno <- convert_genotype(geno_input, 4)

set.seed(2025)
# 50 iterations (id assignment) for all group
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


save(ids,file = "ids_50random among years_5fold_2016_17.RData")

# predict
accuracy=matrix(nrow=length(ids),ncol = 1) #prediction accuracy
pred_y <- list()
file_ids <- list()

for (i in 1:length(ids)) {
    for (j in 1:n.fold) {
      yTrain <- pheno[,1]
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
    file_ids[[i]]<-  data.frame(ID=rownames(pheno),"observed value"=as.numeric(pheno[,1]),GEBV=pred_y_GEBV_all,"predicted value"=pred_all)
    accuracy[i,1]<-cor(as.numeric(pheno[,1]), pred_all, use = "complete.obs")
  }
  
save(file_ids,file = "predicted values for 20167 VIG.RData")
save(accuracy,file = "accuracy for 20167 VIG.RData")

#Genomic prediction: 2016 to predict 2017, 2018
load("/Users/sc3339/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme among years_final/predicted values for 2016 VIG.RData")

accuracy_16vs17=matrix(nrow=50,ncol = 1)
d_2017 <- read.csv("VIG_2017_ar1_BLUE(780).csv",row.names = 1)
for (i in 1:50) {
  identical(rownames(d_2017),file_ids[[i]]$ID)
  accuracy_16vs17[i,1] <- cor(d_2017$VIG, file_ids[[i]]$predicted.value, use = "complete.obs")
}


accuracy_16vs18=matrix(nrow=50,ncol = 1)
d_2018 <- read.csv("VIG_2018_ar1_BLUE(780).csv",row.names = 1)    
for (i in 1:50) {
  identical(rownames(d_2018),file_ids[[i]]$ID)
  accuracy_16vs18[i,1] <- cor(d_2018$VIG, file_ids[[i]]$predicted.value, use = "complete.obs")
}


## Genomic prediction: 2017 to predict 2018

load("/Users/sc3339/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme among years_final/predicted values for 2017 VIG.RData")

accuracy_17vs18=matrix(nrow=50,ncol = 1)
d_2018 <- read.csv("VIG_2018_ar1_BLUE(780).csv",row.names = 1)    
for (i in 1:50) {
  identical(rownames(d_2018),file_ids[[i]]$ID)
  accuracy_17vs18[i,1] <- cor(d_2018$VIG, file_ids[[i]]$predicted.value, use = "complete.obs")
}


## Genomic prediction: 2016 and 2017 to predict 2018

load("/Users/sc3339/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme among years_final/predicted values for 20167 VIG.RData")

accuracy_167vs18=matrix(nrow=50,ncol = 1)
d_2018 <- read.csv("VIG_2018_ar1_BLUE(780).csv",row.names = 1)    
for (i in 1:50) {
  identical(rownames(d_2018),file_ids[[i]]$ID)
  accuracy_167vs18[i,1] <- cor(d_2018$VIG, file_ids[[i]]$predicted.value, use = "complete.obs")
}

### Plot

all <- cbind.data.frame("16vs17"=accuracy_16vs17,"16vs18"=accuracy_16vs18,"17vs18"=accuracy_17vs18,"167vs18"=accuracy_167vs18)

boxplot(all)

all <- rbind.data.frame(accuracy_16vs17,accuracy_16vs18,accuracy_17vs18,accuracy_167vs18)
all$Type <- rep(c("16vs17","16vs18","17vs18","167vs18"),each=50)
all$Type <- factor(all$Type,levels = c("16vs17","16vs18","17vs18","167vs18"))
all$xlab <- c(
  expression(atop("2016", "2017")),
  expression(atop("2016", "2018")),
  expression(atop("2017", "2018")),
  expression(atop("2016&2017", "2018"))
)

box_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728") 

png("boxplot_prediction_amongyears.png",height = 2000,width = 3000,res = 300)
ggplot(all, aes(x = Type, y = V1,fill = Type)) +
  geom_boxplot(size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = box_colors) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(..y.., 3)),  
    color = "red",
    vjust = -0.5,                 
    size = 5
  ) +
  theme_bw() +
  xlab("") +
  ylab("Predictive ability") +
  scale_x_discrete(labels = all$xlab) +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = textGrob(
      expression(atop("Training set", "Testing set")),   # 上下排列
      x = unit(-0.05, "npc"),        
      y = unit(-0.0035, "npc"),         
      just = c("left", "top"),
      gp = gpar(fontsize = 14)
    )
  )
dev.off()

#### SNP marker effect

fit16 <- mixed.solve(y=pheno_16$VIG,Z=t(geno_input),K=NULL,SE=FALSE,return.Hinv = FALSE)

beta16 <- fit16$u


fit17 <- mixed.solve(y=d_2017$VIG,Z=t(geno_input),K=NULL,SE=FALSE,return.Hinv = FALSE)

beta17 <- fit17$u

fit18 <- mixed.solve(y=d_2018$VIG,Z=t(geno_input),K=NULL,SE=FALSE,return.Hinv = FALSE)

beta18 <- fit18$u

fit167 <- mixed.solve(y=pheno$VIG,Z=t(geno_input),K=NULL,SE=FALSE,return.Hinv = FALSE)

beta167 <- fit167$u

marker_cor1 <-cor( beta16,beta18,use="complete")
marker_cor1

marker_cor2 <-cor( beta17,beta18,use="complete")
marker_cor2

marker_cor3 <-cor( beta16,beta17,use="complete")
marker_cor3

marker_cor4 <-cor( beta167,beta18,use="complete")
marker_cor4


library(ggplot2)
library(patchwork)

plot_fun <- function(x, y, xlabel, ylabel){
  
  df <- data.frame(x, y)
  
  r <- cor(x, y, use="complete")
  
  ggplot(df, aes(x, y)) +
    geom_point(
      alpha = 0.4,
      size = 1
    ) +
    geom_smooth(
      method="lm",
      se=FALSE
    ) +
    annotate(
      "text",
      x=Inf,
      y=Inf,
      # label=paste0("r = ", round(r,3)),
      label = paste0("italic(r) == ", round(r, 3)),
      parse = TRUE,
      hjust=1.1,
      vjust=1.5,
      size=6
    ) +
    theme_classic() +
    theme(text = element_text(size = 15))+
    labs(
      x=xlabel,
      y=ylabel
    )
}

p1 <- plot_fun(beta16,beta18,"β2016","β2018")
p2 <- plot_fun(beta17,beta18,"β2017","β2018")
p3 <- plot_fun(beta16,beta17,"β2016","β2017")
p4 <- plot_fun(beta167,beta18,"β2016–17","β2018")

# p1 <- plot_fun(beta16,beta18,"2016","2018")
# p2 <- plot_fun(beta17,beta18,"2017","2018")
# p3 <- plot_fun(beta16,beta17,"2016","2017")
# p4 <- plot_fun(beta167,beta18,"2016–17","2018")

png("/Users/shufen.chen/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/Marker_effects_correlations_across_year2.png",height = 3000,width = 3500,res = 300)
(p3 | p1) /(p2 | p4)
dev.off()


