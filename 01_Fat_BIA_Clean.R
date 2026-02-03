rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(survival)
library(readr)
library(pwr)
library(factoextra)
library(survminer)
library(ggcorrplot)
library(mclust)

source("Codes/00_Basic_Functions.R")

dat_bd <- read_tsv("ukbrapr_data/baseline_dates.tsv")
names(dat_bd)[2:5] <- c(paste0("ins_date_",0:3))
dat_bd_1 <- dat_bd[,c(1:5)]

load("Process_Data/Covar_Data_Inst2.RData")
load("Process_Data/DXA/DXA_BoneMuscle.RData")

dat_abdom <- dat_abdom_u[[1]];abdom_name <- dat_abdom_u[[2]]
abdom_col <- which(stri_detect_fixed(names(dat_abdom),"_20"))
dat_abdom_1 <- dat_abdom[,c(1,abdom_col)]

tar_var <- c("eid",
             "Pancreas_PDFF_fat_fraction_20",
             "Liver_PDFF_fat_fraction_20",
             #"Subcutaneous_fat_volume_20",
             #"Visceral_fat_volume_20",
             #"Visceral_adipose_tissue_volume_VAT_20",
             #"Abdominal_subcutaneous_adipose_tissue_volume_ASAT_20",
             #"Total_trunk_fat_volume_20",
             #"Weight_to_muscle_ratio_20",
             #"Total_abdominal_adipose_tissue_index_20",
             "Muscle_fat_infiltration_20"#,
             #"Posterior_thigh_muscle_fat_infiltration_MFI_left_20",
             #"Posterior_thigh_muscle_fat_infiltration_MFI_right_20",
             #"Anterior_thigh_muscle_fat_infiltration_MFI_right_20",
             #"Anterior_thigh_muscle_fat_infiltration_MFI_left_20",
             #"Abdominal_fat_ratio_20")
)
dat_abdom_1 <- dat_abdom[,which(names(dat_abdom) %in% tar_var)]
dat_abdom_2 <- data.frame(na.omit(dat_abdom_1))

dat_heart <- read_csv("Data_Extract/Data_DL/Heart_MRI_2.csv")
dat_heart_1 <- dat_heart[,c(1,which(stri_detect_fixed(names(dat_heart),"31085")))]
names(dat_heart_1)[2] <- c("area_of_pericardial_fat")

dat.covar %>%
  left_join(dat_abdom_2,by = "eid") %>%
  left_join(dat_heart_1,by = "eid") -> dat.ana

dat.ana <- mutate(dat.ana,height = sqrt(weight/BMI),
                  area_of_pericardial_fat = area_of_pericardial_fat/(height*weight))

dat.ana.1 <- na_remove(data = dat.ana,
                       tar_var = c(tar_var,"sex","area_of_pericardial_fat"))
names(dat.ana.1) <- gsub("_20","",names(dat.ana.1))

save(dat.ana.1,dat.covar,file = "Process_Data/Fat_Dis/Descriptive_Table.RData")

load("Process_Data/Covar_Data_Inst3.RData")

abdom_col_3 <- which(stri_detect_fixed(names(dat_abdom),"_30"))
dat_abdom_3 <- dat_abdom[,c(1,abdom_col_3)]

tar_var_3 <- c("eid","Pancreas_PDFF_fat_fraction_30",
             "Liver_PDFF_fat_fraction_30",
             "Muscle_fat_infiltration_30")

dat_abdom_3_0 <- dat_abdom_3[,which(names(dat_abdom_3) %in% tar_var_3)]
dat_abdom_3_1 <- data.frame(na.omit(dat_abdom_3_0))
dat.covar %>%
  left_join(dat_abdom_3_1,by = "eid") -> dat.ana_3

dat.ana_3_2 <- na_remove(data = dat.ana_3,tar_var = c(tar_var_3[1:3],"sex"))
names(dat.ana_3_2) <- gsub("_30","",names(dat.ana_3_2))

### dat.ana_3_2，包含了脂肪异位和协变量

dat_body <- dat_body_bia_u[[1]];body_name <- dat_body_bia_u[[2]]
body_col_2 <- which(stri_detect_fixed(names(dat_body),"_20"))
body_col_3 <- which(stri_detect_fixed(names(dat_body),"_30"))

dat_body_2 <- dat_body[,c(1,body_col_2)];dat_body_2$instance <- 2
names(dat_body_2) <- gsub("_20","",names(dat_body_2))

dat_body_3 <- dat_body[,c(1,body_col_3)];dat_body_3$instance <- 3
names(dat_body_3) <- gsub("_30","",names(dat_body_3))

identical(names(dat_body_2),names(dat_body_3))

dat_body_fin <- rbind(dat_body_3,dat_body_2)

body_data_um <- dat_body_fin %>%
  mutate(
    Android_Gynoid_ratio = Trunk_fat_mass / ((Leg_fat_mass_left + Leg_fat_mass_right)/2),
    Centrality_index = Trunk_fat_mass / (Arm_fat_mass_left + Arm_fat_mass_right + Leg_fat_mass_left + Leg_fat_mass_right),
    
    Muscle_mass_limbs = (Arm_fat_free_mass_left + Arm_fat_free_mass_right + 
                           Leg_fat_free_mass_left + Leg_fat_free_mass_right),
    Muscle_quality = Muscle_mass_limbs / Weight,
    
    Whole_fat_ratio = Whole_body_fat_free_mass/Whole_body_fat_mass,
    Trunk_fat_ratio = Trunk_fat_free_mass/Trunk_fat_mass,
    Leg_left_fat_ratio = Leg_fat_free_mass_left/Leg_fat_mass_left,
    Leg_right_fat_ratio = Leg_fat_free_mass_right/Leg_fat_mass_right,
    Arm_left_fat_ratio = Arm_fat_free_mass_left/Arm_fat_mass_left,
    Arm_right_fat_ratio = Arm_fat_free_mass_right/Arm_fat_mass_right,    
    Fat_free_mass_percentage = Whole_body_fat_free_mass / Weight * 100,
    Fat_mass_percentage = Whole_body_fat_mass / Weight * 100
  )

muscle_metrics <- c(
  "Whole_body_fat_free_mass", "Trunk_fat_free_mass",
  "Leg_fat_free_mass_left", "Leg_fat_free_mass_right",
  "Arm_fat_free_mass_left", "Arm_fat_free_mass_right"
)

fat_metrics <- c(
  "Body_fat_percentage", "Whole_body_fat_mass", 
  "Trunk_fat_percentage", "Trunk_fat_mass",
  "Leg_fat_percentage_left", "Leg_fat_percentage_right",
  "Leg_fat_mass_left", "Leg_fat_mass_right",
  "Arm_fat_percentage_left", "Arm_fat_percentage_right",
  "Arm_fat_mass_left", "Arm_fat_mass_right","Centrality_index",
  "Whole_fat_ratio","Trunk_fat_ratio","Leg_left_fat_ratio","Leg_right_fat_ratio",
  "Arm_left_fat_ratio","Arm_right_fat_ratio",
  "Fat_mass_percentage",
  "Android_Gynoid_ratio")

dat_body_fin_u <- body_data_um[,c(1,match(c(fat_metrics,muscle_metrics,"instance"),
                                          names(body_data_um)))]

save(dat.ana.1,dat.ana_3_2,dat_body_fin_u,
     file = "Process_Data/Fat_Dis/Fat_Pop.RData")


############   test whether blood variable ################
dat_blood <- read_csv("Data_Extract/Data_DL/Blood_Count.csv")
dat_blood_2 <- dat_blood[,c(1,which(stri_detect_fixed(names(dat_blood),"_i2")))]
dat_blood_2 <- filter(dat_blood_2,!is.na(p30160_i2))

intersect(dat_blood_2$eid,dat.ana.1$eid)