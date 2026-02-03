
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(readr)
library(VIM)
library(caret)
library(splines)
library(glmnet)
library(randomForest)
library(xgboost)
library(pROC)
library(ggpubr)
library(Metrics)
library(plotROC)

source("Codes/00_Basic_Functions.R")
source("Codes_Fat/00_ML_Shared_functions.R")
source("Codes_Fat/00_Fat_BasicFunctions.R")
source("Codes_Fat/00_KOCMI_functions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Fat_Dis/Pro_Impt.RData")

dat_body_fin <- filter(dat_body_fin_u, instance == 2)

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

tar_var_name <- c("age", "sex", "waist_circum", "BMI", tar_var)
dat.ana_2 <- dat.ana.1[, match(c("eid", tar_var_name), names(dat.ana.1))]

dat.ana_2 %>%
  inner_join(pro_int_2_2, by = c("eid")) %>%
  inner_join(dat_body_fin, by = c("eid")) -> dat.ana.u

dat.ana.u <- mutate(dat.ana.u, sex = ifelse(sex == "Male", 1, 0))
dat.ana.u.1 <- fast_impute(dat.ana.u, method = "rf")

na.row <- data.frame(colSums(is.na(dat.ana.u.1)));names(na.row) <- "na_num"

dat.ana.u.1$instance <- NULL
protein_names <- names(dat.ana.u.1)[10:(ncol(dat.ana.u.1) - 27)]
basic_vars <- c("age", "sex", "waist_circum", "BMI")
body_comp_vars <- setdiff(names(dat_body_fin), c("eid", "instance", basic_vars))

###################  for the zero, we used 0.001   ##################
dat.ana.u.1 <- mutate(dat.ana.u.1, sex = as.character(sex))
dat.ana.u.1[,c(6:7)] <- apply(dat.ana.u.1[,c(6:7)],2,function(data){
  log(ifelse(data == 0,0.001,data))
})



# Modified main loop
for (i in c(1:length(tar_var))){
  tar_col_tem <- match(tar_var[i], names(dat.ana.u.1))
  dat_tem <- dat.ana.u.1[, c(tar_col_tem, 2:5, 10:ncol(dat.ana.u.1))]
  names(dat_tem)[1] <- "target"
  load(paste0("Process_Data/Fat_Dis/Continues/", tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$variable_list
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  res_kocmi <- KOCMI.net.parallel(data_processed, k = 3, M = 40)
  save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Fat_Causes/",
                               tar_var[i],"_kocmi.RData"))
}  

# Modified main loop
for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(dat.ana.u.1))
  dat_tem <- dat.ana.u.1[, c(tar_col_tem, 2:5, 10:ncol(dat.ana.u.1))]
  names(dat_tem)[1] <- "target"
  
  load(paste0("Continues/", tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(target_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  for (k_tem in c(3,5,7)){
    res_kocmi <- KOCMI.net.parallel(data_processed, k = k_tem, M = 100)
    save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Fat_Causes/",
                                 tar_var[i],"_",k_tem,"_100_kocmi.RData"))
  }
}  


################################################################################
######################   remove the effects of age and sex  ####################
################################################################################


tar_var_dat <- batch_remove_effects(dat.ana.u.1,outcome_vars= tar_var,
                                    covariate_vars = c("age", "sex"))
names(tar_var_dat) <- gsub("_residual","",names(tar_var_dat))

# Modified main loop
for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(dat.ana.u.1))
  dat_tem <- dat.ana.u.1[, c(tar_col_tem, 2:5, 10:ncol(dat.ana.u.1))]
  dat_tem$target <- tar_var_dat[,i]
  
  load(paste0("Process_Data/Fat_Dis/Continues_Res/", tar_var[i],
               "_Hierarchical_Models.RData"))
  
  selected_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(selected_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  res_kocmi <- KOCMI.net.parallel(data_processed, k = 3, M = 40)
  save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Causes_Res/",tar_var[i],"_kocmi.RData"))
}  

# Modified main loop
for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(dat.ana.u.1))
  dat_tem <- dat.ana.u.1[, c(tar_col_tem, 2:5, 10:ncol(dat.ana.u.1))]
  dat_tem$target <- tar_var_dat[,i]
  
  load(paste0("Process_Data/Fat_Dis/Continues_Res/", 
              tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(target_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  
  for (k_tem in c(3,5,7)){
    res_kocmi <- KOCMI.net.parallel(data_processed, k = k_tem, M = 100)
    save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Causes_Res/",
                                 tar_var[i],"_",k_tem,"100_kocmi.RData"))
  }
}  

################################################################################
###########  remove the effects of age and sex on protein  ####################
################################################################################

tar_var_all <- batch_remove_effects(dat.ana.u.1,
                                    outcome_vars= names(dat.ana.u.1)[-c(1:5)],
                                    covariate_vars = c("age", "sex"),
                                    model_type = "lm")
names(tar_var_all) <- gsub("_residual","",names(tar_var_all))
tar_var_all <- data.frame(apply(tar_var_all,2,scale))

tar_var_nm <- batch_remove_effects(dat.ana.u.1,
                                    outcome_vars= names(dat.ana.u.1)[-c(1:5)],
                                    covariate_vars = c("age", "sex"),
                                    nonliner_vars = "age",
                                    model_type = "nlm")
names(tar_var_nm) <- gsub("_residual","",names(tar_var_nm))
tar_var_nm <- data.frame(apply(tar_var_nm,2,scale))

save(tar_var_all,tar_var,tar_var_nm,file = "Fat_KOCMI.RData")

# Modified main loop
for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(tar_var_all))
  dat_tem <- tar_var_all[, c(tar_col_tem, 5:ncol(tar_var_all))]
  names(dat_tem)[1] <- "target"
  
  load(paste0("Process_Data/Fat_Dis/Continues/", 
              tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(target_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  
  for (k_tem in c(3,5,7)){
    res_kocmi <- KOCMI.net.parallel(data_processed, k = k_tem, M = 100)
    save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Fat_Causes_Nor_Res/",
                                 tar_var[i],"_",k_tem,"100_kocmi.RData"))
  }
}  


for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(tar_var_all))
  dat_tem <- tar_var_all[, c(tar_col_tem, 5:ncol(tar_var_all))]
  names(dat_tem)[1] <- "target"
  
  load(paste0("Process_Data/Fat_Dis/Continues/", 
              tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(target_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  
  dat_model <- dat_tem[,match(c("target",target_proteins),names(dat_tem))]
  data_processed <- preprocess_data(data = dat_model)
  
  for (k_tem in c(3,5,7)){
    res_kocmi <- KOCMI.net(data_processed, k = k_tem, M = 100)
    save(res_kocmi,file = paste0("Process_Data/Fat_Dis/Fat_Causes_Nor_Res/",
                                 tar_var[i],"_",k_tem,"100_kocmi.RData"))
  }
}  


################################################################################
###########  Pearson correlation coefficient for directions ####################
################################################################################

pro_mer <- c()

for (i in c(1:length(tar_var))){
  print(i)
  tar_col_tem <- match(tar_var[i], names(tar_var_all))
  dat_tem <- tar_var_all[, c(tar_col_tem, 5:ncol(tar_var_all))]
  names(dat_tem)[1] <- "target"
  
  load(paste0("Process_Data/Fat_Dis/Continues/", 
              tar_var[i], "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$selected_proteins
  
  nonlinear_results <- analyze_nonlinear_effects(target_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  pro_mer <- c(pro_mer,target_proteins)
}  

dat_model <- tar_var_all[,c(1:4,match(unique(pro_mer),names(tar_var_all)))]
cor_results <- psych::corr.test(dat_model, method = "pearson", adjust = "none")
cor_results <- as.data.frame(cor_results$r) %>%
  tibble::rownames_to_column(var = "var1") %>%
  pivot_longer(cols = -var1, names_to = "var2", values_to = "correlation")

save(cor_results,file = "Process_Data/Fat_Dis/Fat_Causes_Nor_Res/Correlation.RData")

