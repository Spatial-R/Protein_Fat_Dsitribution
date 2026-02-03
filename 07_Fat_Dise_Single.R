
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(survival)
library(readr)
library(factoextra)
library(survminer)
library(cowplot)
library(ggcorrplot)
library(mclust)
source("Codes/00_Basic_Functions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Disease.RData")
load("Process_Data/Covar_Data_Inst2.RData")
dat_death <- read_tsv("ukbrapr_data/death.tsv")
dat_death <- mutate(dat_death,eid = as.numeric(eid))

dat.covar <- dat.ana.1

dat_body_fin <- filter(dat_body_fin_u, instance == 2)

tar_var <- c("Muscle_fat_infiltration","Liver_PDFF_fat_fraction",
             "Pancreas_PDFF_fat_fraction","area_of_pericardial_fat")

dir_tar <- "Process_Data/Cases_Full"
file_dis <- list.files(path = dir_tar,pattern = "csv",full.name = T)[c(18,21,74)]

res_fin <- data.frame()

model_vars <- list(
  model_2_var = model_2_var, model_3_var = model_3_var[1:6]
)

for (i in c(1:length(file_dis))) {
  
  print(i)
  dat_dem <-  read.csv(file_dis[i],stringsAsFactors = F)
  names_tem <- gsub("/","",gsub(".csv","",gsub(dir_tar,"",file_dis[i])))
  tar_col <- match(c("eid","df","bin","bin_prev"),names(dat_dem))
  
  dat.covar %>%
    left_join(dat_dem[,tar_col],by = "eid") %>%
    left_join(dat_death[,c(2,6)],by = "eid") -> dat.ana

  dat.ana.1 <- mutate(dat.ana,
                      age = as.numeric(substr(ins_date2,1,4)) - birth_year,
                      df = ifelse(bin == 0 & !is.na(date_of_death),as.character(date_of_death),as.character(df)),
                      t_time = as.numeric(as.Date(df) - as.Date(ins_date2)))
  
  names(dat.ana.1)[match(c("df","bin","bin_prev"),names(dat.ana.1))] <- 
    c("out_df","out_bin","out_bin_prev")
  
  dat.ana.2 <- filter(dat.ana.1,out_bin_prev == 0)
  #dat.ana.2 <- filter(dat.ana.1,!(out_bin == 1 & (t_time < 2*365)))
  names(dat.ana.2)[match(c("out_bin"),names(dat.ana.2))] <- c("status")
  
  dat.ana.3 <- na_remove(data = dat.ana.2,
                         tar_var = c(tar_var,model_2_var,model_3_var[2:6]))
  tar_cols <- match(tar_var,names(dat.ana.3))
  dat.ana.3[,tar_cols] <- apply(dat.ana.3[,tar_cols],2,scale)
  
  dat_res_multi <- regression_fit(data = dat.ana.3,tar_var = tar_var[1],
                              ori_full = T,rem_var = NULL,
                              add_var = tar_var[-1],ref_level = NULL,
                              exp_var_list = list(model_2_var,model_3_var[2:6]),
                              out_full =TRUE)
  dat_res_multi <- mutate(dat_res_multi,type = "multi")
  dat_res_single <- data.frame()
  
  for (k in tar_var){
    dat_res_tem <- regression_fit(data = dat.ana.3,tar_var = k,
                                ori_full = T,rem_var = NULL,
                                add_var = NULL,ref_level = NULL,
                                exp_var_list = list(model_2_var,model_3_var[2:6]),
                                out_full =TRUE)   
    dat_res_tem <- mutate(dat_res_tem,type = "single")
    dat_res_single <- rbind(dat_res_single,dat_res_tem)
  }
  dat_res_2 <- rbind(dat_res_single,dat_res_multi)
  dat_res_2 <- mutate(dat_res_2,id = names_tem);
  res_fin <- rbind(res_fin,dat_res_2)
}


res_fin_u1 <- mutate(res_fin,p_ind_u = ifelse(p_ind < 0.05,1,0))
res_fin_u1 <- filter(res_fin_u1,var %in% tar_var)

res_fin_u3 <- filter(res_fin_u1, p_ind_u == 1) 
write.csv(res_fin_u1,file = "Process_Data/Fat_Dis/Fat_Dis_Sig_1.csv",row.names = F)
write.csv(res_fin_u1,file = "/Users/zhangbing/Documents/Tmp/Fat_Dis_1.csv",row.names = F)
# write.csv(res_fin_u2[,c(6,8,2,4)],file = "Process_Data/Fat_Dis/Fat_Dis_All.csv",row.names = F)

label_mn <- c("Pancreas PDFF","Liver PDFF","Muscle FI","Pericardial fat area")
target_dis <- c("diabetes_type_2","coronary_heart","stroke_all")
target_dis_u <- c("Type 2 Diabetes","Coronary Heart Disease","Stroke")

res_fin_u2 <- mutate(res_fin_u1,HR = exp(coef),
                     Low_95 = exp(coef - 1.96*se), 
                     High_95 = exp(coef + 1.96*se),
                     type = factor(type,levels = c("single","multi"),
                                   labels = c("Individually adjusted models",
                                              "Mutually adjusted models")),
                     var = factor(var,levels = tar_var, labels = label_mn),
                     model_type = factor(model_type,levels = c(0:3),
                                         labels = c(paste0("Model ",0:3))),
                     id = factor(id,levels = target_dis,labels = target_dis_u))


fig_1 <- ggplot(data = res_fin_u2) +
  geom_pointrange(aes(x = var, y = HR,ymin = Low_95, ymax = High_95,
                      color = model_type,
                      group = model_type),
                  size = 0.2,position = position_dodge2(width = 0.5)) +
  facet_grid(type~id,scale = "free_y") + 
  scale_color_manual(values = cols[c(3,6,8)],
                     guide = guide_legend(title = "Models: ")) +
  geom_hline(yintercept = 1,linetype = 2,color = "black") +
  xlab("") + ylab("Hazard ratio (95%CI)") +
  theme_bw(base_size = 12,base_family = "serif") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(filename = "Figures/Fat/Fat_Diseases.tiff",fig_1,
       width = 22,height = 14,units = "cm",dpi = 300)
