library(dplyr)
library(tidyverse)
library(ggrepel)
load("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/accuracy_group for scheme2.RData")

scheme1 <- accuracy_cycle

load("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/accuracy_group for scheme3.RData")

scheme2 <- accuracy_cycle


Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")

combined_data_scheme1 <- data.frame()  
for (i in 1:length(scheme1)) {
  result_scheme1 <- do.call(rbind, scheme1[[i]])
  combined_data_scheme1 <- rbind.data.frame(combined_data_scheme1, result_scheme1)
}

combined_data_scheme2 <- data.frame()  
for (i in 1:length(scheme2)) {
  result_scheme2 <- do.call(rbind, scheme2[[i]])
  combined_data_scheme2 <- rbind.data.frame(combined_data_scheme2, result_scheme2)
}

combined_data_scheme1$group <- rep(c("OTTM","SIBR","CASIA","EURO","Check"),each=6)

combined_data_scheme1$scheme <- rep(paste0(seq(100,0,by=-20),"%"),times=250)

combined_data_scheme2$group <- rep(c("OTTM","SIBR","CASIA","EURO","Check"),each=6)

combined_data_scheme2$scheme <- rep(paste0(seq(100,0,by=-20),"%"),times=250)

combined_data_melt_scheme1 <- reshape::melt(combined_data_scheme1,id.vars=c("group","scheme"))
combined_data_melt_scheme2 <- reshape::melt(combined_data_scheme2,id.vars=c("group","scheme"))
#write.table(combined_data,file = "combined_data for scheme4_OTTM.txt",row.names = F)
#combined_data_order <- combined_data[order(combined_data$scheme,decreasing = T),]

combined_data_melt_scheme1$group <- factor(combined_data_melt_scheme1$group,levels=c("OTTM","SIBR","CASIA","EURO","Check"))
combined_data_melt_scheme2$group <- factor(combined_data_melt_scheme2$group,levels=c("OTTM","SIBR","CASIA","EURO","Check"))

combined_data_melt_all <- rbind.data.frame(combined_data_melt_scheme1,combined_data_melt_scheme2)
combined_data_melt_all$Senario <- rep(c("Scheme1","Scheme2"),each=4500)

# combined_data_summary <- combined_data_melt_all %>%
#   group_by(across(c(1:3,5))) %>%
#   summarise(
#     mean = mean(.data[["value"]], na.rm = TRUE),
#     sd = sd(.data[["value"]], na.rm = TRUE),
#     .groups = "drop"
#   )

# write.csv(combined_data_summary ,"~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/Scheme1&2_mean_PV_SD.csv",row.names = F)

### line add error bar
combined_data_melt_all$scheme <- factor(combined_data_melt_all$scheme,levels = paste0(seq(0,100,by=20),"%"))

label_data <- combined_data_melt_all %>%
  # filter(scheme %in% c("0%", "100%")) %>%
  filter(scheme %in% c( "100%")) %>%
  group_by(Senario,scheme, variable, group) %>%
  summarise(mean_value = mean(value), .groups = "drop")

png("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/Scheme1&2_line.png",height = 3000,width = 3500,res = 300)
ggplot(combined_data_melt_all, aes(x = scheme, y = value, 
                                   group = interaction(Senario,variable, group),
                                  colour =Senario )) +
  scale_colour_manual(values = c("#999999","#E69F00"))+
  # scale_shape_manual(values = c(3,5))+
  stat_summary(fun = mean, geom = "line" ) +
  stat_summary(fun = mean, geom = "point", size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15) +
  
  # stat_summary(fun = mean, geom = "line", 
  #              position = position_dodge(width = 0.3)) +
  # 
  # stat_summary(fun = mean, geom = "point", 
  #              position = position_dodge(width = 0.3), size = 2) +
  
  # stat_summary(fun.data = mean_se, geom = "errorbar", 
  #              width = 0.15, 
  #              position = position_dodge(width = 0.3)) +
  
  geom_text_repel(data = label_data,
                  aes(y = mean_value, label = round(mean_value, 3)),
                  position = position_dodge(width = 0.3),
                  vjust = -0.9,
                  hjust=1,
                  fontface = "bold",
                  size = 4,
                  show.legend = FALSE)+
  facet_grid(
    factor(group, levels = c("OTTM","SIBR","CASIA","EURO","Check")) ~
      factor(variable, levels = c("GH","HGT","VIG"))
  ) +
  theme_bw() +
  xlab("") +
  ylab("Average predictive ability") +
  # labs(fill = "Group", color = "Group") +
  theme(text = element_text(size = 15), legend.position = "right")

dev.off()

#### Scheme 1 Kinship mean
load("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/file_group for scheme2.RData")
file_scheme1 <- file_cycle
Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
res_scheme1 <- data.frame()
percentage_train <- paste0("train", "_",seq(50,0,by=-10))
for (i in 1:50) {
  for (g in 1:length(Group)) {
    for (j in 1:length(percentage_train)) {
      kin_mean= mean(file_scheme1[[i]][[g]][["GH"]][[j]]$k_mean)
      res_scheme1 <- rbind(res_scheme1, data.frame(
        cycle = i,
        group = Group[g],
        train=percentage_train[j],
        kin_mean = kin_mean))
      
    }
  }
}
res_scheme1$scheme <- rep(paste0(seq(100,0,by=-20),"%"),times=250)

res_scheme1$group[res_scheme1$group=="55H94|Genoa|Hybriforce3400|Vernal"]<- "Check"

kin_ability_scheme1 <- cbind.data.frame(res_scheme1,combined_data_scheme1)
kin_ability_scheme12 <- kin_ability_scheme1[c(2,4,5:8)]
kin_ability_scheme13 <- reshape::melt(kin_ability_scheme12,id.vars=c("group","scheme","kin_mean"))
kin_ability_scheme13$group[kin_ability_scheme13$group=="55H94|Genoa|Hybriforce3400|Vernal"]<- "Check"

kin_ability_scheme_summary <- kin_ability_scheme13 %>%
  group_by(across(c(1:2,4))) %>%
  summarise(
    mean = mean(.data[["value"]], na.rm = TRUE),
    mean_kin = mean(.data[["kin_mean"]], na.rm = TRUE),
    sd = sd(.data[["value"]], na.rm = TRUE),
    sd_kin = sd(.data[["kin_mean"]], na.rm = TRUE),
    .groups = "drop"
  )

#### Scheme 2 Kinship mean
load("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/file_group for scheme3.RData")
file_scheme2 <- file_cycle
Group <- c("OTTM","SIBR","CASIA","EURO","55H94|Genoa|Hybriforce3400|Vernal")
res_scheme2 <- data.frame()
percentage_target <- paste0("target", "_",seq(50,0,by=-10))
for (i in 1:50) {
  for (g in 1:length(Group)) {
    for (j in 1:length(percentage_target)) {
      kin_mean= mean(file_scheme2[[i]][[g]][["GH"]][[j]]$k_mean)
      res_scheme2 <- rbind(res_scheme2, data.frame(
        cycle = i,
        group = Group[g],
        train=percentage_target[j],
        kin_mean = kin_mean))
      
    }
  }
}
res_scheme2$scheme <- rep(paste0(seq(100,0,by=-20),"%"),times=250)

res_scheme2$group[res_scheme2$group=="55H94|Genoa|Hybriforce3400|Vernal"]<- "Check"

kin_ability_scheme21 <- cbind.data.frame(res_scheme2,combined_data_scheme2)
kin_ability_scheme22 <- kin_ability_scheme21[c(2,4,5:8)]
kin_ability_scheme23 <- reshape::melt(kin_ability_scheme22,id.vars=c("group","scheme","kin_mean"))
kin_ability_scheme23$group[kin_ability_scheme23$group=="55H94|Genoa|Hybriforce3400|Vernal"]<- "Check"

kin_ability_scheme2_summary <- kin_ability_scheme23 %>%
  group_by(across(c(1:2,4))) %>%
  summarise(
    mean = mean(.data[["value"]], na.rm = TRUE),
    mean_kin = mean(.data[["kin_mean"]], na.rm = TRUE),
    sd = sd(.data[["value"]], na.rm = TRUE),
    sd_kin = sd(.data[["kin_mean"]], na.rm = TRUE),
    .groups = "drop"
  )

## all kin
all_kin_ability <- rbind.data.frame(kin_ability_scheme13,kin_ability_scheme23)
all_kin_ability$Senario <- rep(c("Scheme1","Scheme2"),each=4500)
# write.csv(all_kin_ability,"all_kin_ability.csv",row.names = F)

plot_data_kin <- all_kin_ability %>%
  group_by(Senario, variable, group, scheme) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    mean_kin = mean(kin_mean, na.rm = TRUE),
    se_value = sd(value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )



# plot_data_kin %>%
#   group_by(group) %>%
#   summarise(
#     xmin = min(mean_kin_round),
#     xmax = max(mean_kin_round)
#   )

# plot_data_kin$mean_kin <- factor(round(plot_data_kin$mean_kin,3), levels = sort(unique(round(plot_data_kin$mean_kin,3)),decreasing = F))
plot_data_kin$scheme <- factor(plot_data_kin$scheme,levels = paste0(seq(0,100,by=20),"%"))

plot_data_kin$mean_kin_round <- round(plot_data_kin$mean_kin, 4)
summary(plot_data_kin$mean_kin_round)


png("~/Desktop/Alfalfa_Revision1_Due_2026-05-04/Genomic prediction/results2/Scheme1&2_kinship&ability2.png",height = 3000,width = 4500,res = 300)
ggplot(plot_data_kin,aes(x = mean_kin_round,y = mean_value,colour = Senario,group = Senario)) +
  # scale_colour_manual(values = c("#999999","#E69F00"))+
  scale_colour_manual(values = c("#56B4E9", "#CC79A7"))+
  geom_line(linewidth=0.1) +
  # scale_x_continuous(breaks = seq(-0.015,-0.005,0.001))+
  # scale_x_continuous(limits = c(-0.017,-0.002))+
  # scale_x_continuous(
  #   breaks = sort(unique(plot_data_kin$mean_kin_round))
  # ) +
  geom_point(size = 1) +
  # geom_text(aes(label = scheme),
  #           vjust = -0.8,
  #           size = 3,show.legend = FALSE) +
  # geom_text_repel(aes(label = scheme),
  #                 size = 2,
  #                 show.legend = FALSE)+
  geom_text_repel(
    aes(label = scheme),
    size = 5,
    show.legend = FALSE,
    max.overlaps = Inf
  )+

  geom_errorbar(aes(ymin = mean_value - se_value,
                    ymax = mean_value + se_value)) +

  facet_grid(
    factor(variable, levels = c("GH","HGT","VIG")) ~
      factor(group, levels = c("OTTM","SIBR","CASIA","EURO","Check")), scales = "free_x") +
  scale_x_continuous(
    labels = function(x) sprintf("%.4f", x)
  )+

  theme_bw()+
  
  # ggtitle("The scatter plot of mean kinship and average predictive ability \n by adding additional samples (0%–100%, in 20% increments) to training set")+
  xlab("Genetic relatedness between the test and training sets") +
  # xlab("Mean kinship of the test set to the training set at each proportion ((0%–100%, in 20% increments))") +
  ylab("Average predictive ability") +
  theme(
    text = element_text(size = 20),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

dev.off()

