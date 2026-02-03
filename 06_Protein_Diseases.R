
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
library(magrittr)
library(doSNOW)

source("Codes/00_Basic_Functions.R")
load("Process_Data/Disease.RData")
load("Process_Data/Covar_Data.RData")

dat_death <- read_tsv("ukbrapr_data/death.tsv")
dat_death <- mutate(dat_death,eid = as.numeric(eid))

dat_protein_0 <- read_csv("Data_Extract/Data_DL/proteomic.csv")
names(dat_protein_0) <- gsub("olink_instance_0.","",names(dat_protein_0))

all_icd_name <- c(paste0("D",stri_pad_left(0:99,width = 2,"0")),
                  paste0("E",stri_pad_left(0:99,width = 2,"0")),
                  paste0("F",stri_pad_left(0:99,width = 2,"0")),
                  paste0("G",stri_pad_left(0:99,width = 2,"0")),
                  paste0("I",stri_pad_left(0:99,width = 2,"0")),
                  paste0("J",stri_pad_left(0:99,width = 2,"0")),
                  paste0("K",stri_pad_left(0:99,width = 2,"0")),
                  paste0("M",stri_pad_left(0:99,width = 2,"0")),
                  paste0("N",stri_pad_left(0:30,width = 2,"0")))
dis_icd <- mutate(dis_icd,icd = stri_trim_both(icd,pattern = "\\P{Wspace}"))
dis_icd_0 <- filter(dis_icd,icd %in% all_icd_name)
dis_sum <- data.frame(table(dis_icd_0$icd))
dis_sum_0 <- filter(dis_sum,Freq > 1000)

loop_dat <- expand.grid(icd = as.character(dis_sum_0$Var1),
            pro = names(dat_protein_0)[-1])


mcors <- 60
cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
set.seed(998468235L,kind = "L'Ecuyer")
mcor.set <- list(preschedule = FALSE,set.seed = TRUE)

global.result <- foreach(i = 1:nrow(loop_dat), .combine=rbind,
                         .packages = c("dplyr","magrittr","survival"),.errorhandling = "remove",
                         .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
                           
                           icd_name <- loop_dat[i,1];pro_name <- loop_dat[i,2]
                           dat_dem <- filter(dis_icd_0,icd == as.character(icd_name))
                           names(dat_dem)[2] <- "dis_date"
                           
                           dat_protein_0 %>%
                             left_join(dat_dem[,1:2],by = "eid") %>%
                             left_join(dat_death[,c(2,6)],by = "eid") %>%
                             left_join(dat.covar[,-c(26:82)],by = "eid") -> dat.ana
                           
                           dat.ana.1 <- mutate(dat.ana,
                                               age = as.numeric(substr(ins_date0,1,4)) - birth_year,
                                               bin = ifelse(!is.na(dis_date),1,0),
                                               df = dplyr::if_else(is.na(dis_date) & !is.na(date_of_death),
                                                                   date_of_death,dis_date),
                                               df = dplyr::if_else(!is.na(dis_date),dis_date,as.Date("2024-10-30")),
                                               t_time = as.numeric(as.Date(df) - as.Date(ins_date0))/365.25)
                           dat.ana.2 <- filter(dat.ana.1,t_time > 0)
                           
                           names(dat.ana.2)[match(c("bin"),names(dat.ana.2))] <- c("status")
                           names(dat.ana.2)[match(as.character(pro_name),names(dat.ana.2))] <- c("target")
                           
                           case_num <- length(which(dat.ana.2$status == 1))
                           
                           dat.ana.3 <- na_remove(data = dat.ana.2,
                                                  tar_var = c("target",model_2_var,model_3_var[1:6]))
                           
                           dat_res_2 <- regression_fit(data = dat.ana.3,tar_var = "target",
                                                           ori_full = T,rem_var = NULL,
                                                           add_var = NULL,ref_level = NULL,
                                                           exp_var_list = list(model_2_var,model_3_var[1:6]),
                                                           out_full =TRUE)
                          write.csv(dat_res_2,file = paste0(pro_name,"*",icd_name,"_reg.csv"))
                          dat_res_3 <- filter(dat_res_2,var == "target")
                          dat_res_3 <- mutate(dat_res_3,pop_num = nrow(dat.ana.3),
                                              case_num = case_num,icd = icd_name,protein = pro_name)
                          return(dat_res_3)     
                            
  }
save(global.result,file = "Pro_Icd_Reg.RData")            
stopCluster(cl)
registerDoSEQ()


load("/Users/zhangbing/Downloads/Pro_Icd_Reg.RData")
global.result <- mutate(global.result,icd = as.character(icd))
dis_name_1 <- mutate(dis_name_1,des_1 = stri_trim_both(des_1,pattern = "\\P{Wspace}"))
global.result.1 <- merge(global.result,dis_name_1,by.x = "icd",by.y = "des_1")
global.result.2 <- filter(global.result.1,model_type == 2 & p_ind < 0.05)
View(dat_mn <- filter(global.result.2,protein %in% c("lrrn1")))
write.csv(dat_mn[,c(2,10)],file = "mx.csv",row.names = F)






res_fin <- data.frame(); dat.covar$egfr <- NULL

for (i in c(1:nrow(dis_sum_0))) {
  print(i)
  dat_dem <- filter(dis_icd,icd == dis_sum_0[i,1])
  names(dat_dem)[2] <- "dis_date"

  dat_protein_0 %>%
    left_join(dat_dem[,1:2],by = "eid") %>%
    left_join(dat_death[,c(2,6)],by = "eid") %>%
    left_join(dat.covar[,-c(26:82)],by = "eid") -> dat.ana

  dat.ana.1 <- mutate(dat.ana,
                      age = as.numeric(substr(ins_date0,1,4)) - birth_year,
                      bin = ifelse(!is.na(dis_date),1,0),
                      df = dplyr::if_else(is.na(dis_date) & !is.na(date_of_death),
                                          date_of_death,dis_date),
                      df = dplyr::if_else(!is.na(dis_date),dis_date,as.Date("2024-10-30")),
                      t_time = as.numeric(as.Date(df) - as.Date(ins_date0))/365.25)
  dat.ana.2 <- filter(dat.ana.1,t_time > 0)

  names(dat.ana.2)[match(c("bin"),names(dat.ana.2))] <- c("status")

  dat.ana.3 <- na_remove(data = dat.ana.2,
                         tar_var = c(model_2_var,model_3_var[1:6]))

  for (j in c(2:ncol(dat_protein_0))){
    dat.ana.4 <- dat.ana.3
    tar_col <- match(names(dat_protein_0)[j],names(dat.ana.4))
    if(!is.na(tar_col)){
      names(dat.ana.4)[match(names(dat_protein_0)[j],names(dat.ana.4))] <- "target"
      dat.ana.5 <- filter(dat.ana.4,!is.na(target))

      dat_res_2 <- regression_fit(data = dat.ana.5,tar_var = "target",
                                  ori_full = T,rem_var = NULL,
                                  add_var = NULL,ref_level = NULL,
                                  exp_var_list = list(model_2_var,model_3_var[1:6]),
                                  out_full =TRUE)
      dat_res_2 <- mutate(dat_res_2,num = nrow(dat.ana.5),id = dis_sum[i,1],
                          protein = names(dat_protein_0)[j])
      res_fin <- rbind(res_fin,dat_res_2)
    }
  }
}

res_fin_u1 <- filter(res_fin,var %in% "target" & model_type == 2)
res_fin_u1 <- mutate(res_fin_u1,p_ind_u = ifelse(p_ind < 0.05,1,0))
res_fin_u2 <- merge(res_fin_u1,dis_name_1,by.x = 'id',by.y = "des_1")
res_fin_u3 <- filter(res_fin_u2, p_ind_u == 1)
write.csv(res_fin_u3[,c(1,6,8,2,4)],file = "Results_sig.csv",row.names = F)
write.csv(res_fin_u2[,c(6,8,2,4)],file = "Results.csv",row.names = F)


################################################################################
#############################  Subset of diseases ##############################
################################################################################

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

loop_tar <- c()

for ( k in c(1:length(tar_var))){
  print(k)
  load(file = paste0("Process_Data/Fat_Dis/Continues/", tar_var[k], "_Hierarchical_Models.RData"))
  target_pro <- all_models_results[[1]]$variable_list
  loop_tar <- c(loop_tar,target_pro)
}  
  

  dir_tar <- "Process_Data/Cases_Full"
  file_dis <- list.files(path = dir_tar,pattern = "csv",full.name = T)[c(18,21,74)]
  
  loop_dat <- expand.grid(icd = as.character(file_dis),
                          pro = unique(loop_tar))
  
  mcors <- 2
  cl <- makeCluster(mcors,type = "SOCK"); registerDoSNOW(cl)
  set.seed(998468235L,kind = "L'Ecuyer")
  mcor.set <- list(preschedule = FALSE,set.seed = TRUE)
  
  global.result <- foreach(i = 1:nrow(loop_dat), .combine=rbind,
                           .packages = c("dplyr","magrittr","survival"),.errorhandling = "remove",
                           .inorder = FALSE, .options.multicore = mcor.set) %dopar%  {
                             
                             icd_name <- as.character(loop_dat[i,1]);
                             pro_name <- as.character(loop_dat[i,2]);
                             
                             dat_dem <-  read.csv(icd_name,stringsAsFactors = F)
                             names_tem <- gsub("/","",gsub(".csv","",gsub(dir_tar,"",icd_name)))
                             tar_col <- match(c("eid","df","bin","bin_prev"),names(dat_dem))
                             
                             dat_protein_0 %>%
                               left_join(dat_dem[,tar_col],by = "eid") %>%
                               left_join(dat_death[,c(2,6)],by = "eid") %>%
                               left_join(dat.covar[,-c(26:82)],by = "eid") -> dat.ana
                             
                             dat.ana.1 <- mutate(dat.ana,
                                                 age = as.numeric(substr(ins_date0,1,4)) - birth_year,
                                                 df = ifelse(bin == 0 & !is.na(date_of_death),as.character(date_of_death),as.character(df)),
                                                 t_time = as.numeric(as.Date(df) - as.Date(ins_date0)))
                             
                             names(dat.ana.1)[match(c("df","bin","bin_prev"),names(dat.ana.1))] <- 
                               c("out_df","out_bin","out_bin_prev")
                             
                             dat.ana.2 <- filter(dat.ana.1,out_bin_prev == 0)
                             
                             names(dat.ana.2)[match(c("out_bin"),names(dat.ana.2))] <- c("status")
                             names(dat.ana.2)[match(as.character(pro_name),names(dat.ana.2))] <- c("target")
                             
                             case_num <- length(which(dat.ana.2$status == 1))
                             
                             dat.ana.3 <- na_remove(data = dat.ana.2,
                                                    tar_var = c("target",model_2_var,model_3_var[1:6]))
                             
                             dat_res_2 <- regression_fit(data = dat.ana.3,tar_var = "target",
                                                         ori_full = T,rem_var = NULL,
                                                         add_var = NULL,ref_level = NULL,
                                                         exp_var_list = list(model_2_var,model_3_var[1:6]),
                                                         out_full =TRUE)
                             write.csv(dat_res_2,file = paste0("Process_Data/Pro_Dis/",pro_name,"*",names_tem,"_reg.csv"))
                             dat_res_3 <- filter(dat_res_2,var == "target")
                             dat_res_3 <- mutate(dat_res_3,pop_num = nrow(dat.ana.3),
                                                 case_num = case_num,icd = names_tem,protein = pro_name)
                             return(dat_res_3)     
                           }
  save(global.result,file = "Process_Data/Pro_Diseases_CardMet.RData")            
  stopCluster(cl)
  registerDoSEQ()

  
  target_dis <- c("diabetes_type_2","coronary_heart","stroke_all")
  target_dis_u <- c("Type 2 Diabetes","Coronary Heart Disease","Stroke")
  
  for ( k in c(1:length(tar_var))){
    print(k)
    load(file = paste0("Process_Data/Fat_Dis/Continues/", tar_var[k], "_Hierarchical_Models.RData"))
    target_pro <- all_models_results[[1]]$variable_importance$average$Protein
    res_tem <- filter(global.result,protein %in% target_pro & model_type == 2)
    
    res_tem_1 <- mutate(res_tem,HR = exp(coef),
                         Low_95 = exp(coef - 1.96*se), 
                         High_95 = exp(coef + 1.96*se),
                         p_ind_u = ifelse(p_ind < 0.05,1,0),
                        protein = factor(protein,levels = target_pro),
                         icd = factor(icd,levels = target_dis,labels = target_dis_u))
    fig_tem <- ggplot(data = res_tem_1) +
      geom_pointrange(aes(x = protein, y = HR,ymin = Low_95, ymax = High_95,
                          color = factor(p_ind_u)),
                      size = 0.1,position = position_dodge2(width = 0.5)) +
      facet_grid(~icd,scale = "free_y") + 
      coord_flip()+
      scale_color_manual(values = cols[c(6,8)],
                         guide = "none") +
      geom_hline(yintercept = 1,linetype = 2,color = "black") +
      xlab("") + ylab("Hazard ratio (95%CI)") +
      theme_bw(base_size = 12,base_family = "serif") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "bottom")
    
    ggsave(filename = paste0("Figures/Fat/Pro_Dis/",tar_var[k],".tiff"),fig_tem,
           width = 17,height = 14,units = "cm",dpi = 300)
  }  
