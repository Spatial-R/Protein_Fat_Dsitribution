rm(list = ls())

library(dplyr)
library(stringi)
library(tidyr)
library(readr)
library(VIM)

re_run <- FALSE

############# 肌肉中脂肪组织浸润只测量了一次，肝脏和胰腺的PDFF测了有2次 ######

source("Codes/00_Basic_Functions.R")

load("Process_Data/Fat_Dis/Fat_Pop.RData")

protein_2 <- read_csv("Data_Extract/Data_DL/proteomic_2.csv")
names(protein_2) <- gsub("olink_instance_2.","",names(protein_2))
protein_3 <- read_csv("Data_Extract/Data_DL/proteomic_3.csv")
names(protein_3) <- gsub("olink_instance_3.","",names(protein_3))

pro_name <- read.csv("Data/Protein/UKB_PPP.csv",header = T)
pro_name <- mutate(pro_name,meaning = tolower(meaning),
                   meaning = gsub("-","_",meaning))

int_eid_1 <- intersect(protein_2$eid,dat.ana.1$eid);
int_eid_2 <- intersect(protein_3$eid,dat.ana_3_2$eid);

if(isTRUE(re_run)){
    pro_int_2_0 <- filter(protein_2,eid %in% int_eid_1)
    na.row.prot_2 <- data.frame(colSums(is.na(pro_int_2_0)));names(na.row.prot_2) <- "na_num"
    rm_dat_2 <- filter(na.row.prot_2,na_num > 90)
    pro_int_2_1 <-  pro_int_2_0[,-c(match(row.names(rm_dat_2),names(pro_int_2_0)))]
    pro_int_2_2 <- fast_impute(pro_int_2_1,method = "rf") ### imputing missing values  ####
    
    pro_int_3_0 <- filter(protein_3,eid %in% int_eid_2)
    na.row.prot_3 <- data.frame(colSums(is.na(pro_int_3_0)));names(na.row.prot_3) <- "na_num"
    rm_dat_3 <- filter(na.row.prot_3,na_num > 90)
    pro_int_3_1 <-  pro_int_3_0[,-c(match(row.names(rm_dat_3),names(pro_int_3_0)))]
    pro_int_3_2 <- fast_impute(pro_int_3_1,method = "rf")   ### imputing missing values ####
    save(pro_int_3_2,pro_int_2_2,file = "Process_Data/Fat_Dis/Pro_Impt.RData")
} else {
    load("Process_Data/Fat_Dis/Pro_Impt.RData")
}

identical(names(pro_int_3_2),names(pro_int_2_2))
tar_pro_name <- intersect(names(pro_int_3_2),names(pro_int_2_2))

pro_int_3_3 <- pro_int_3_2[,match(tar_pro_name,names(pro_int_3_2))]
pro_int_2_3 <- pro_int_2_2[,match(tar_pro_name,names(pro_int_2_2))]
pro_int_3_3$instance <- 3; pro_int_2_3$instance <- 2;
pro_int_fin <- rbind(pro_int_2_3,pro_int_3_3)

ukb_pro_sel <- filter(pro_name,meaning %in% names(pro_int_2_2))

setdiff(names(pro_int_2_2),pro_name$meaning)

save(pro_int_fin,file = "Process_Data/Fat_Dis/Protein_Cleaned.RData")
write.csv(ukb_pro_sel,file = "Results/Fat/Table S2.csv",row.names = F)
