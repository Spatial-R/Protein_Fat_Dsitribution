
rm(list = ls())

library(dplyr)
library(stringi)
library(tidyr)
library(tableone)

source("Codes/00_Basic_Functions.R")

#load("Process_Data/Fat_Dis/Descriptive_Table.RData")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Fat_Dis/Pro_Impt.RData")

res_vars <- c("age","sex","rural_category","townsend_index",
              "edu_type","waist_circum", "BMI_cut","smoke_status",
               "drink_status","phy_act",'sleep_cut')

dat.ana.2 <- filter(dat.ana.1,! Liver_PDFF_fat_fraction == 0)
inter_id <- intersect(dat.ana.2$eid,pro_int_2_2$eid)

dat.ana.full <- dat.ana.2[,match(c("eid",res_vars),names(dat.ana.2))]
dat.ana.tar <- filter(dat.ana.full,eid %in% inter_id)
dat.ana.full <- mutate(dat.ana.full,type = "full")
dat.ana.tar <- mutate(dat.ana.tar,type = "target")
dat.ana.fin <- rbind(dat.ana.tar,dat.ana.full)
dat.ana.fin$eid <- NULL

dat.ana.fin <- mutate(dat.ana.fin,
                      BMI_cut = factor(BMI_cut,levels = c("Underweight","Normal weight",
                                                    "Overweight","Obesity")),
                      sleep_cut = factor(sleep_cut,levels = c("Less","Normal","More")),
                      drink_status = factor(drink_status,levels = c("Never","Previous","Current")),
                      smoke_status = factor(smoke_status,levels = c("Never","Previous","Current")),
                      phy_act  = factor(phy_act, levels = c("Low","Moderate","High")))
                      
                      
tab1 <- CreateTableOne(vars = res_vars,data = dat.ana.fin,
                       factorVars = res_vars[c(2:3,5,8:11)],
                       strata ="type")
tab1mat <- print(tab1,quote = FALSE,noSpace = TRUE,printToggle = TRUE)

write.csv(tab1mat,file = "Results/Fat/Table S3.csv",row.names = T)


