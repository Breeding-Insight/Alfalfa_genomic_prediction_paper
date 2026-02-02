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

set.seed(2025)
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
set.seed(123)
# Initialize lists to store results
cycles <- 50
file_percen <- list()
pheno_train_file <- list()
geno_train_file <- list()
accuracy <- matrix(NA, nrow=5, ncol=length(pheno_all))
rownames(accuracy) <- paste0("train", "_",seq(50,10,by=-10))
colnames(accuracy) <- names(pheno_all)

file_pheno <- list()
pheno_train_pheno <- list()
geno_train_pheno <- list()

file_cycle <- list()
pheno_train_cycle <- list()
geno_train_cycle <- list()
accuracy_cycle <- list()

file_group <- list()
pheno_train_group <- list()
geno_train_group <- list()
accuracy_group <- list()

  for (r in 1:cycles) {
    for (i in 1:length(Group)) {
    #scheme4
   
      test <- as.matrix(sample(1:nrow(pheno_group[[i]]),round(0.5*nrow(pheno_group[[i]]))))
      train_50 <- setdiff(1:nrow(pheno_group[[i]]),test) #50%
      train_40 <- sample(train_50, round(0.4 * nrow(pheno_group[[i]]))) #40%
      train_30 <- sample(train_40, round(0.3 * nrow(pheno_group[[i]]))) #30%
      train_20 <- sample(train_30, round(0.2 * nrow(pheno_group[[i]]))) #20%
      train_10 <- sample(train_20, round(0.1 * nrow(pheno_group[[i]]))) #10%
      
      percentage <- list(train_50, train_40, train_30, train_20, train_10)
      names(percentage) <- paste0("train", "_",seq(50,10,by=-10))
      
      for (p in 1:length(pheno_all)) {
      for (j in 1:length(percentage)) {
        pheno_train_all <- rbind(pheno_all[!rownames(pheno_all) %in% rownames(pheno_group[[i]]), ],
                                 pheno_group[[i]][percentage[[j]], ])
        geno_train_all <- rbind(geno_all[!rownames(geno_all) %in% colnames(geno_group[[i]]), ],
                                t(geno_group[[i]])[percentage[[j]], ])
        print(identical(rownames(pheno_train_all),rownames(geno_train_all)))
        pheno_test <- pheno_group[[i]][test, ]
        geno_test <- t(geno_group[[i]])[test, ]
        
        y <- pheno_train_all[, p]
        y_model <- mixed.solve(y, Z=geno_train_all, K=NULL, SE=FALSE, return.Hinv=FALSE)
        e <- as.matrix(y_model$u)
        pred_y_GEBV <- geno_test %*% e
        pred_y <- pred_y_GEBV[, 1] + as.numeric(y_model$beta)
        pred_y <- data.frame(pred_y)
        pred_y <- pred_y[rownames(pheno_test), ]
        pred_y_GEBV_all <- pred_y - as.numeric(y_model$beta)
        y_valid <- pheno_test[, p]
        
        file_percen[[j]] <- data.frame(ID=rownames(pheno_test), "observed value"=y_valid,GEBV=pred_y_GEBV_all,
                                       "predicted value"=pred_y)
        names(file_percen)[j] <- names(percentage)[j]
        pheno_train_file[[j]] <- pheno_train_all
        names(pheno_train_file)[j] <- names(percentage)[j]
        geno_train_file[[j]] <- geno_train_all
        names(geno_train_file)[j] <- names(percentage)[j]
        accuracy[j, p] <- cor(pred_y, y_valid, use="complete")
        
      }
      
      file_pheno[[p]] <- file_percen
      names(file_pheno)[p] <- names(pheno_all)[p]
      pheno_train_pheno[[p]] <- pheno_train_file
      names(pheno_train_pheno)[p] <- names(pheno_all)[p]
      geno_train_pheno[[p]] <- geno_train_file
      names(geno_train_pheno)[p] <- names(pheno_all)[p]
    }
      file_group[[i]] <- file_pheno
      names(file_group)[i] <- Group[i]
      pheno_train_group[[i]] <- pheno_train_pheno
      names(pheno_train_group)[i] <- Group[i]
      geno_train_group[[i]] <- geno_train_pheno
      names(geno_train_group)[i] <- Group[i]
      accuracy_group[[i]] <- accuracy 
      names(accuracy_group)[i] <- Group[i]
      
  }
  
    file_cycle[[r]] <- file_group
    pheno_train_cycle[[r]] <- pheno_train_group
    geno_train_cycle[[r]] <- geno_train_group
    accuracy_cycle[[r]] <- accuracy_group
  
}
save(file_cycle,file = "file_group for scheme4_update.RData")
save(pheno_train_cycle,file = "pheno_train_group for scheme4_update.RData")
#save(geno_train_group,file = "geno_train_group for scheme4_update.RData")
save(accuracy_cycle,file = "accuracy_group for scheme4_update.RData")

load("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme4/accuracy_group for scheme4_update.RData")
Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")

combined_data <- data.frame()  
for (i in 1:length(accuracy_cycle)) {
  result <- do.call(rbind, accuracy_cycle[[i]])
  combined_data <- rbind.data.frame(combined_data, result)
}
combined_data$group <- rep(c("OTTM","SIBR","CASIA","EURO","Check"),each=5)
# combined_data$scheme <- rep(paste0("train", "_",seq(50,10,by=-10)),times=5)
combined_data$scheme <- rep(paste0(seq(50,10,by=-10),"%"),times=5)
combined_data_melt <- reshape::melt(combined_data,id.vars=c("group","scheme"))
#write.table(combined_data,file = "combined_data for scheme4_OTTM.txt",row.names = F)
#combined_data_order <- combined_data[order(combined_data$scheme,decreasing = T),]

combined_data_melt$group <- factor(combined_data_melt$group,levels=c("OTTM","SIBR","CASIA","EURO","Check"))


library(dplyr)

combined_data_summary <- combined_data_melt %>%
  group_by(across(1:3)) %>%
  summarise(
    mean = mean(.data[["value"]], na.rm = TRUE),
    sd = sd(.data[["value"]], na.rm = TRUE),
    .groups = "drop"
  )

write.csv(combined_data_summary ,"Scheme4_mean_PV_SD.csv",row.names = F)

##### plus Scheme1
load("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme1/predicted values for scheme1.rda")
#accuracy for each group
accuracy <- list()
for (i in 1:length(Group)) {
  accuracy[[i]] <- list()
  for (j in 1:length(file_all)) {
    group_value <- file_all[[j]][grep(pattern = Group[i],file_all[[j]]$ID),]
    accuracy[[i]][[j]] <- cor(group_value$observed.value,group_value$predicted.value, use = "complete.obs")
    names(accuracy[[i]])[[j]]<- names(file_all)[j]
  }
  names(accuracy)[[i]]<- Group[i]
}

d_unlist <- data.frame(value=unlist(accuracy))
d_unlist$variable <- stringr::str_split_fixed(rownames(d_unlist),pattern = "\\.",n=2)[,2]
d_unlist$group <- stringr::str_split_fixed(rownames(d_unlist),pattern = "\\.",n=2)[,1]
d_unlist$group[which(d_unlist$group==Group[5])] <- "Check"
d_unlist$scheme <- "0%(Scheme1)"
d_unlist$scheme <- "0%"
data_melt <- d_unlist[,c("group","scheme","variable","value")]

# all_melt <- rbind.data.frame(data_melt,combined_data_mean)
all_melt <- rbind.data.frame(data_melt,combined_data_melt)
all_melt$group <- factor(all_melt$group,levels=c("OTTM","SIBR","CASIA","EURO","Check"))
# all_melt$scheme[all_melt$scheme=="50%"] <- "50%"


## add error bar
library(ggplot2)
library(ggrepel)

label_data <- all_melt %>%
  filter(scheme %in% c("0%(Scheme1)", "50%")) %>%
  group_by(scheme, variable, group) %>%
  summarise(mean_value = mean(value), .groups = "drop")

png("~/Desktop/Alfalfa_GWAS/20250313/Genomic prediction/Results/scheme4/Scheme4_line_update_4.png",height = 2000,width = 5700,res = 300)
ggplot(all_melt, aes(x = scheme, y = value, 
                     group = interaction(variable, group),
                     fill = group, color = group)) +
  stat_summary(fun = mean, geom = "line", 
               position = position_dodge(width = 0.3)) +

  stat_summary(fun = mean, geom = "point", 
               position = position_dodge(width = 0.3), size = 2) +

  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, 
               position = position_dodge(width = 0.3)) +
  geom_text_repel(data = label_data,
            aes(y = mean_value, label = round(mean_value, 3)),
            position = position_dodge(width = 0.3),
            vjust = -0.9,
            hjust=1,
            fontface = "bold",
            size = 4,
            show.legend = FALSE)+
  facet_wrap(~factor(variable, levels = c("GH","HGT","VIG"))) +
  theme_bw() +
  xlab("") +
  ylab("Value") +
  labs(fill = "Group", color = "Group") +
  theme(text = element_text(size = 15), legend.position = "right")
dev.off()

