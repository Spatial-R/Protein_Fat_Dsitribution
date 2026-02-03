
library(TwoSampleMR)
library(dplyr)
library(ggplot2)
#usethis::edit_r_environ()

gcst_ids <- paste0("ebi-a-",c(
  "GCST90016666", "GCST90016667", "GCST90016668", "GCST90016669", "GCST90016670",
  "GCST90016671", "GCST90016672", "GCST90016673", "GCST90016674", "GCST90016675", 
  "GCST90016676"))[c(10,8)]

#### ebi-a-GCST90016673   Percent liver fat

res_clump <- read.csv("Data/Protein/Clumped.csv")

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

dat_res <- data.frame()  ; # length(tar_var)

for (i in c(1:2)){
  
  load(paste0("Process_Data/Fat_Dis/Continues_Res/", tar_var[i],
              "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$variable_list
  target_proteins <- toupper(target_proteins)
  
   res_clump_1 <- filter(res_clump,Symbol %in% target_proteins)
   loop_var <- unique(res_clump_1$Symbol)
   res_list <- lapply(loop_var, function(data){
     print(match(data,loop_var))
     dat_tem <- filter(res_clump_1,Symbol == data)
     repeat{
       try({sle <- extract_outcome_data(snps = dat_tem$SNP,outcomes = gcst_ids[i])})
       if(exists("sle")){break}
          Sys.sleep(0.5)
     }
     ana_dat_tem <- harmonise_data(exposure_dat = dat_tem,outcome_dat = sle,action = 2)
     mr_res <- mr(ana_dat_tem);
     if(nrow(mr_res) == 0){
       mr_res <- data.frame() 
     }else{
       mr_res <- mutate(mr_res,id = data) 
     }
     return(mr_res)
   })
   mr_res_all <- bind_rows(res_list)
   mr_res_all <- mutate(mr_res_all,var = tar_var[i])
   dat_res <- rbind(dat_res,mr_res_all)
}
dat_res_u <- filter(dat_res,method == "Inverse variance weighted" & pval < 0.05)

View(filter(dat_res,method == "Inverse variance weighted" & pval < 0.05))


dat_res_sen <- data.frame()  ; # length(tar_var)

for (i in c(1:2)){
  
  load(paste0("Process_Data/Fat_Dis/Continues/", tar_var[i],
              "_Hierarchical_Models.RData"))
  
  target_proteins <- all_models_results[[1]]$variable_list
  target_proteins <- toupper(target_proteins)
  
  res_clump_1 <- filter(res_clump,Symbol %in% target_proteins)
  loop_var <- unique(res_clump_1$Symbol)
  res_list <- lapply(loop_var, function(data){
    print(match(data,loop_var))
    dat_tem <- filter(res_clump_1,Symbol == data)
    repeat{
      try({sle <- extract_outcome_data(snps = dat_tem$SNP,outcomes = gcst_ids[i])})
      if(exists("sle")){break}
      Sys.sleep(0.5)
    }
    ana_dat_tem <- harmonise_data(exposure_dat = dat_tem,outcome_dat = sle,action = 2)
    mr_res <- mr(ana_dat_tem);
    if(nrow(mr_res) == 0){
      mr_res <- data.frame() 
    }else{
      mr_res <- mutate(mr_res,id = data) 
    }
    return(mr_res)
  })
  mr_res_all <- bind_rows(res_list)
  mr_res_all <- mutate(mr_res_all,var = tar_var[i])
  dat_res_sen <- rbind(dat_res_sen,mr_res_all)
}
save(dat_res_sen,dat_res,file = "Process_Data/Fat_Dis/Protein_MR.RData")

dat_res_sen_u <- filter(dat_res_sen,method == "Inverse variance weighted" & pval < 0.05)
View(dat_res_sen_u)
