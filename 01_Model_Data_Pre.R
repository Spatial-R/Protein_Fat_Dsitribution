
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
source("Codes_Fat/00_ML_Classified_functions.R")
source("Codes_Fat/00_Fat_BasicFunctions.R")

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
dat.ana.u.1 <- filter(dat.ana.u.1,! Liver_PDFF_fat_fraction == 0)
dat.ana.u.1[,c(6:7)] <- apply(dat.ana.u.1[,c(6:7)],2,function(data){
  log(ifelse(data == 0,0.001,data))
})

###############################  Classified  ##################################

cardmus_thsd <- dat.ana.u.1 %>%
  group_by(sex) %>%
  summarise(
    thr_car = quantile(area_of_pericardial_fat, 0.75, na.rm = TRUE),
    thr_muc = quantile(Muscle_fat_infiltration, 0.75, na.rm = TRUE),.groups = 'drop')

# Define gender-specific clinical thresholds
clinical_thresholds <- list(
  # Male thresholds
  Male = list(
    Pancreas_PDFF = 10.4,  Liver_PDFF = 5.0,     
    Muscle = as.numeric(cardmus_thsd[2,3]), 
    area_of_peri = as.numeric(cardmus_thsd[2,2])),
  Female = list(
    Pancreas_PDFF = 10.4,   Liver_PDFF = 5.0,     
    Muscle = as.numeric(cardmus_thsd[1,3]), 
    area_of_peri = as.numeric(cardmus_thsd[1,2])))

dat_classfied <- create_binary_outcomes(data = dat.ana.u.1,
                                      thresholds = clinical_thresholds)

save(dat.ana.u.1,dat_classfied,file = "Process_Data/Fat_Dis/Model_Data.RData")
