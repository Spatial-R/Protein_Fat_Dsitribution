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

setwd("/public/home/eyyangchongzhe/hpc1-data/UKB/")

load("Server.RData")

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
                           write.csv(dat_res_2,file = paste0("Full_Res/",pro_name,"*",icd_name,"_reg.csv"))
                           dat_res_3 <- filter(dat_res_2,var == "target")
                           dat_res_3 <- mutate(dat_res_3,pop_num = nrow(dat.ana.3),
                                               case_num = case_num,icd = icd_name,protein = pro_name)
                           return(dat_res_3)     
                           
                         }
save(global.result,file = "Pro_Icd_Reg.RData")            
stopCluster(cl)
registerDoSEQ()
