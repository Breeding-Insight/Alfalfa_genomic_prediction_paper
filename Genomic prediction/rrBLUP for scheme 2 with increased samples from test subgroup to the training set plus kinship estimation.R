library(rrBLUP)
library(tidyverse)


set.seed(2025)
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
identical(rownames(geno_all),rownames(pheno_all))

K <- A.mat(geno_all) 

#genomic prediction

# Initialize lists to store results
cycles <- 50
file_percen <- list()
pheno_train_file <- list()
geno_train_file <- list()
accuracy <- matrix(NA, nrow=6, ncol=length(pheno_all))
rownames(accuracy) <- paste0("target", "_",seq(100,0,by=-20))
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
    #Scheme3
    test <- as.matrix(sample(1:nrow(pheno_group[[i]]),round(0.5*nrow(pheno_group[[i]]))))
    target_50 <- setdiff(1:nrow(pheno_group[[i]]),test) #100% from the remaining 50% test subgroup
    target_40 <- sample(target_50, round(0.4 * nrow(pheno_group[[i]]))) #80%
    target_30 <- sample(target_40, round(0.3 * nrow(pheno_group[[i]]))) #60%
    target_20 <- sample(target_30, round(0.2 * nrow(pheno_group[[i]]))) #40%
    target_10 <- sample(target_20, round(0.1 * nrow(pheno_group[[i]]))) #20%
    target_0 <- sample(target_10, round(0 * nrow(pheno_group[[i]]))) #0%
    
    percentage_target <- list(target_50, target_40, target_30, target_20, target_10,target_0)
    names(percentage_target) <- paste0("target", "_",seq(100,0,by=-20))
    
    size_train <- nrow(pheno_all[!rownames(pheno_all) %in% rownames(pheno_group[[i]]), ])
    train<- sample(1:size_train,size_train-length(target_50))
    
    for (p in 1:length(pheno_all)) {
      
      for (j in 1:length(percentage_target)) {
        pheno_train_all <- rbind(pheno_all[!rownames(pheno_all) %in% rownames(pheno_group[[i]]), ][train,],
                                 pheno_group[[i]][percentage_target[[j]], ])
        geno_train_all <- rbind(geno_all[!rownames(geno_all) %in% colnames(geno_group[[i]]), ][train,],
                                t(geno_group[[i]])[percentage_target[[j]], ])
        print(identical(rownames(pheno_train_all),rownames(geno_train_all)))
        pheno_test <- pheno_group[[i]][test, ]
        geno_test <- t(geno_group[[i]])[test, ]
        # identical(rownames(geno_test),rownames(pheno_test))
        #kinship
        geno_kin <- rbind(geno_train_all, geno_test)
        geno_kin2 <- K[which(rownames(K)%in% rownames(geno_kin)),which(colnames(K)%in% rownames(geno_kin))]
        geno_kin3 <- geno_kin2[rownames(geno_kin), rownames(geno_kin)]
        identical(rownames(geno_kin3),rownames(geno_kin))
        identical(colnames(geno_kin3),rownames(geno_kin))
        
        n_train <- nrow(geno_train_all)
        K_test_train <- geno_kin3[(n_train+1):nrow(geno_kin3), 1:n_train]
        identical(rownames(geno_test),rownames(K_test_train))
        K_avg <- data.frame(K_mean=rowMeans(K_test_train))
        # identical(rownames(K_avg),rownames(pheno_test))
        K_avg <- K_avg[rownames(pheno_test), ]
        
        #####
        y <- pheno_train_all[, p]
        y_model <- mixed.solve(y, Z=geno_train_all, K=NULL, SE=FALSE, return.Hinv=FALSE)
        e <- as.matrix(y_model$u)
        pred_y_GEBV <- geno_test %*% e
        pred_y <- pred_y_GEBV[, 1] + as.numeric(y_model$beta)
        pred_y <- data.frame(pred_y)
        # identical(rownames(pred_y),rownames(pheno_test))
        pred_y <- pred_y[rownames(pheno_test), ]
        pred_y_GEBV_all <- pred_y - as.numeric(y_model$beta)
        y_valid <- pheno_test[, p]
        
        file_percen[[j]] <- data.frame(ID=rownames(pheno_test), "observed value"=y_valid,GEBV=pred_y_GEBV_all,
                                       "predicted value"=pred_y, k_mean=K_avg)
        names(file_percen)[j] <- names(percentage_target)[j]
        pheno_train_file[[j]] <- pheno_train_all
        names(pheno_train_file)[j] <- names(percentage_target)[j]
        geno_train_file[[j]] <- geno_train_all
        names(geno_train_file)[j] <- names(percentage_target)[j]
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
