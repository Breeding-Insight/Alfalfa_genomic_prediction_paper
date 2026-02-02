setwd("~/Desktop/Alfalfa_GWAS/20250313/Pheno")
asreml.license.activate()
detach("package:asreml", unload=TRUE)
###############################################################################
##############Mixed Model analysis of Alfalfa with Repeated Measurements ######
##############Prepared by Bhoja Basnet, Breeding Insight  ###################
###############################################################################
rm(list=ls())

######################Loading libraries
library(tidyverse)
#library(dplyr)
library(asreml)
library(ggplot2)
asreml.options(maxit=20)

############################Reading data
alf<-read.csv("VIG",header=TRUE)
alf <- alf[grep(pattern = "2016|2017",alf$year),]
#alf <- subset(alf,alf$year=="2018")
head(alf)
dim(alf)
str(alf)
unique(alf$date)
######## Defining factors###############
alf$date <- factor(alf$date)
alf$year <- factor(alf$year)
alf$cut <- factor(alf$cut)
alf$rep <- factor(alf$rep)
alf$blk <- factor(alf$blk)
alf$row <- factor(alf$row)
alf$rng <- factor(alf$rng)
alf$plant_no <-factor(alf$plant_no)
alf$plot <- factor(alf$plot)
alf$name <- factor(alf$name)
alf$plant_id <-factor(alf$plant_id)
alf$plant_name <- factor(alf$plant_name)
#alf$check<- factor(alf$check)
alf$id_by_plot<-factor(alf$id_by_plot)
alf$id_by_plant_aug<-factor(alf$id_by_plant_aug)


###### Checking the structure of the dataset 
str(alf)

###############################################################################
###############BLUE with FIXED EFFECT plant_id Repeated Measures Model#########
###############################################################################
alf<-alf[order(alf$plant_id, alf$cut),] #Make sure it is in order as in residual terms
alf<-alf[order(alf$cut, alf$plant_id),] #Make sure it is in order

asreml.options(pworkspace = "5gb")

modAR1_hom_fixed_plants <- asreml(fixed= value ~ cut + plant_id,
                                  random = ~ at(cut):rng + at(cut):row,  # AR1 for rows and columns with homogeneous variance
                              residual = ~ ar1(cut):id(plant_id),# AR1 for cuts (with heterogeneous variance for cuts)
                              data = alf
                              )

modAR1_hom_fixed_plants <-update.asreml(modAR1_hom_fixed_plants)

###########Fixed effect model 2

modAR2_hom_fixed_plants <- asreml(fixed= value ~ cut + plant_id,
                                  random = ~ at(cut):rng:row,  # AR1 for rows and columns with homogeneous variance
                                  residual = ~ ar1(cut):id(plant_id),# AR1 for cuts (with heterogeneous variance for cuts)
                                  data = alf
                                  )

modAR2_hom_fixed_plants <-update.asreml(modAR2_hom_fixed_plants)

###########Fixed effect model 3
modAR3_hom_fixed_plants <- asreml(fixed= value ~ cut + plant_id,
                                  random = ~ at(cut):ar1(rng):ar1(row),  # AR1 for rows and columns with homogeneous variance
                                  residual = ~ ar1(cut):id(plant_id),# AR1 for cuts (with heterogeneous variance for cuts)
                                  data = alf
                                  )


modAR3_hom_fixed_plants <-update.asreml(modAR3_hom_fixed_plants)

############## Model evaluation- comparing three models

loglik1<-summary(modAR1_hom_fixed_plants)$loglik
loglik2<-summary(modAR2_hom_fixed_plants)$loglik
loglik3<-summary(modAR3_hom_fixed_plants)$loglik
loglik1
loglik2
loglik3
#BIC
summary(modAR1_hom_fixed_plants)$bic[1]
summary(modAR2_hom_fixed_plants)$bic[1]
summary(modAR3_hom_fixed_plants)$bic[1]
#AIC
summary(modAR1_hom_fixed_plants)$aic[1]
summary(modAR2_hom_fixed_plants)$aic[1]
summary(modAR3_hom_fixed_plants)$aic[1]

lrt_1_vs_2 <- 2 * (loglik2 - loglik1)
p_value_1_vs_2 <- pchisq(lrt_1_vs_2, df = 1, lower.tail = FALSE)
cat("Model 1 vs Model 2: Chi-square =", lrt_1_vs_2, " P-value =", p_value_1_vs_2, "\n")

lrt_2_vs_3 <- 2 * (loglik3 - loglik2)
p_value_2_vs_3 <- pchisq(lrt_2_vs_3, df = 1, lower.tail = FALSE)
cat("Model 2 vs Model 3: Chi-square =", lrt_2_vs_3, " P-value =", p_value_2_vs_3, "\n")

####modAR2_hom_fixed_plants is better than modAR1_hom_fixed_plants 

plot(modAR3_hom_fixed_plants)
summary(modAR3_hom_fixed_plants)$varcomp
wald.asreml(modAR3_hom_fixed_plants)
pred_BLUE_repeat<-predict(modAR3_hom_fixed_plants, classify="cut:plant_id")
head(pred_BLUE_repeat)
pred_BLUE_repeat_combined<-predict(modAR3_hom_fixed_plants, classify="plant_id")
head(pred_BLUE_repeat_combined)

#suppressWarnings(write.csv(pred_BLUE_repeat,file="VIG_2016_2017_ar1_BLUE of plants_repeated measurement combined by_cut_plant.csv", row.names = FALSE))
suppressWarnings(write.csv(pred_BLUE_repeat_combined,file="VIG_2016_2017_ar1_BLUE of plants_repeated measurement combined by_plant.csv", row.names = FALSE))

####################Heritability- combined analysis#############################
#repeatability: BLUP based on AIC/BIC selection for BLUE (across year)

#GH VIG
modAR3_hom_random_plants <- asreml(fixed= value ~ cut ,
                                  random = ~ at(cut):ar1(rng):ar1(row)+ plant_id,  
                                  residual = ~ ar1(cut):id(plant_id),
                                  data = alf
)

modAR3_hom_random_plants <-update.asreml(modAR3_hom_random_plants)
summary(modAR3_hom_random_plants)$varcomp

vc <- summary(modAR3_hom_random_plants)$varcomp
vc_with_index <- vc %>% # Add index numbers
  mutate(Index = row_number())
print(vc_with_index)

G_Var <-vc$component[13] #plant_id

pred <- suppressWarnings(asreml::predict.asreml(modAR3_hom_random_plants, classify = "plant_id", sed = TRUE, trace = 0))
pred$pvals

vdBLUP.mat <- pred$sed^2 # variance of difference of predicted "plant_id" values
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])#average var of diff = sq of std err diff
vdBLUP.avg
h2_cullis <- 1 - vdBLUP.avg /(2*G_Var) #H2 Estimate proposed by Cullis et al. 2006 
h2_cullis #(choose)



#HGT
modAR2_hom_random_plants <- asreml(fixed= value ~ cut ,
                                  random = ~ at(cut):rng:row + plant_id,  # AR1 for rows and columns with homogeneous variance
                                  residual = ~ ar1(cut):id(plant_id),# AR1 for cuts (with heterogeneous variance for cuts)
                                  data = alf
)

summary(modAR2_hom_random_plants)$varcomp

vc <- summary(modAR2_hom_random_plants)$varcomp
vc_with_index <- vc %>% # Add index numbers
  mutate(Index = row_number())
print(vc_with_index)

G_Var <-vc$component[6] #plant_id

pred <- suppressWarnings(asreml::predict.asreml(modAR2_hom_random_plants, classify = "plant_id", sed = TRUE, trace = 0))
pred$pvals

vdBLUP.mat <- pred$sed^2 # variance of difference of predicted "plant_id" values
vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag = FALSE)])#average var of diff = sq of std err diff
vdBLUP.avg
h2_cullis <- 1 - vdBLUP.avg /(2*G_Var) #H2 Estimate proposed by Cullis et al. 2006 
h2_cullis #(choose)

