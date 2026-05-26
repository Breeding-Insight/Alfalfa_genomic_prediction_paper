library(ggpubr)
library(rstatix)
library(GWASpoly)
setwd("~/Desktop/Alfalfa_GWAS/20250313/PCA&Kinship")


#ALL
set.seed(2015)
data <- read.GWASpoly(ploidy=4, pheno.file="Pheno(780).csv", geno.file="geno_promis_bias_od_maf_filer for all groups_GWAS with maf >0.03.csv",
                      format="numeric", n.traits=3, delim=",")

data.loco <- set.K(data,LOCO=TRUE,n.core=2)
N <- 780 #Population size
params <- set.params(geno.freq = 1 - 5/N)

data.loco.scan <- GWASpoly(data=data.loco,models=c("additive"),
                           traits="GH",params=params,n.core=10)
data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)


####
p <- LD.plot(data2)
p + xlim(0,30)+ylim(0,0.2)

p[["data"]][p[["data"]]$r2==0.1,]

#OTTM
set.seed(2015)

data_OTTM <- read.GWASpoly(ploidy=4, pheno.file="Pheno_OTTM(178).csv", geno.file="geno_promis_bias_od_maf_filer for OTTM_GWAS with maf >0.03.csv",
                      format="numeric", n.traits=3, delim=",")

data.loco_OTTM <- set.K(data_OTTM ,LOCO=TRUE,n.core=2)
N_OTTM <- 178 #Population size
params_OTTM <- set.params(geno.freq = 1 - 5/N_OTTM)

data.loco.scan_OTTM <- GWASpoly(data=data.loco_OTTM,models=c("additive"),
                           traits="GH",params=params,n.core=2)
data2_OTTM <- set.threshold(data.loco.scan_OTTM,method="M.eff",level=0.05)

p_OTTM <- LD.plot(data2_OTTM)
p_OTTM + xlim(0,30) 


#SIBR
set.seed(2015)

data_SIBR <- read.GWASpoly(ploidy=4, pheno.file="Pheno_SIBR(176).csv", geno.file="geno_promis_bias_od_maf_filer for SIBR_GWAS with maf >0.03.csv",
                          format="numeric", n.traits=3, delim=",")


data.loco_SIBR<- set.K(data_SIBR ,LOCO=TRUE,n.core=2)
N_SIBR <- 176 #Population size
params_SIBR <- set.params(geno.freq = 1 - 5/N_SIBR)

data.loco.scan_SIBR <- GWASpoly(data=data.loco_SIBR,models=c("additive"),
                                traits="GH",params=params,n.core=2)
data2_SIBR <- set.threshold(data.loco.scan_SIBR,method="M.eff",level=0.05)

p_SIBR <- LD.plot(data2_SIBR)

#CASIA
set.seed(2015)

data_CASIA <- read.GWASpoly(ploidy=4, pheno.file="Pheno_CASIA(155).csv", geno.file="geno_promis_bias_od_maf_filer for CASIA_GWAS with maf >0.03.csv",
                           format="numeric", n.traits=3, delim=",")


data.loco_CASIA <- set.K(data_CASIA ,LOCO=TRUE,n.core=2)
N_CASIA <- 155 #Population size
params_CASIA <- set.params(geno.freq = 1 - 5/N_CASIA)

data.loco.scan_CASIA <- GWASpoly(data=data.loco_CASIA,models=c("additive"),
                                traits="GH",params=params,n.core=2)
data2_CASIA <- set.threshold(data.loco.scan_CASIA,method="M.eff",level=0.05)

p_CASIA <- LD.plot(data2_CASIA)

#EURO
set.seed(2015)
data_EURO <- read.GWASpoly(ploidy=4, pheno.file="Pheno_EURO(180).csv", geno.file="geno_promis_bias_od_maf_filer for EURO_GWAS with maf >0.03.csv",
                            format="numeric", n.traits=3, delim=",")


data.loco_EURO <- set.K(data_EURO ,LOCO=TRUE,n.core=2)
N_EURO <- 180 #Population size
params_EURO <- set.params(geno.freq = 1 - 5/N_EURO)

data.loco.scan_EURO <- GWASpoly(data=data.loco_EURO,models=c("additive"),
                                 traits="GH",params=params,n.core=2)
data2_EURO <- set.threshold(data.loco.scan_EURO,method="M.eff",level=0.05)

p_EURO <- LD.plot(data2_EURO)


#Check
set.seed(2015)
data_Check <- read.GWASpoly(ploidy=4, pheno.file="Pheno_Check(91).csv", geno.file="geno_promis_bias_od_maf_filer for Check_GWAS with maf >0.03.csv",
                           format="numeric", n.traits=3, delim=",")


data.loco_Check <- set.K(data_Check ,LOCO=TRUE,n.core=2)
N_EURO <- 91 #Population size
params_Check <- set.params(geno.freq = 1 - 5/N_EURO)

data.loco.scan_Check <- GWASpoly(data=data.loco_Check,models=c("additive"),
                                traits="GH",params=params,n.core=2)
data2_Check <- set.threshold(data.loco.scan_Check,method="M.eff",level=0.05)

p_Check <- LD.plot(data2_Check)

data_all <- rbind.data.frame(p$data,p_OTTM$data,p_SIBR$data,p_CASIA$data,p_EURO$data,p_Check$data)
#write.csv(data_all,"LD.csv",row.names = F)
data_all$Group <- rep(c("All","OTTM","SIBR","CASIA","EURO","Check"),each=500)
data_all$Group <- factor(data_all$Group,levels =c("All","OTTM","SIBR","CASIA","EURO","Check") )

png("~/Desktop/Alfalfa_GWAS/20250313/PCA&Kinship/LD.png",width=3500, height=3000,res = 300)
ggplot(data = data_all, aes(x = data_all$d, y = data_all$r2,col=Group)) + ylab(expression(r^2)) + 
  geom_hline(yintercept = 0.1, linetype="dashed",col="grey",size=1 )+
  # annotate("text",x=0,y=0.1,label="0.025",col="red",size=5)+
  # ylim(c(0,1))+
  scale_color_manual(values=c("#addd8e","#fdae6b","#9ebcda","#8856a7","#e34a33","#1c9099"))+
  theme_bw() + geom_line(size=1,alpha=0.9)+ ylim(0,0.15)+xlim(0,20)+xlab("Distance (Mb)")+theme(text = element_text(size=20))
dev.off()

library(dplyr)

data_all %>%
  group_by(Group) %>%
  summarise(
    d_at_r2_0.1 = approx(r2, d, xout = 0.1)$y
  )
