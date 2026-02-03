
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(igraph)
library(ggraph)

source("Codes_Fat/00_Causes_Network_functions.R")
tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

tar_var_abb <- c("Pancreatic PDFF","Liver PDFF","Muscle FI","PFA")

load("Process_Data/Fat_Dis/Fat_Causes_Nor_Res/Correlation.RData")
names(cor_results)[1:2] <- c("regulator","target")

cor_results_1 <- cor_results %>%
  mutate(across(1:2, ~ recode(., "Pancreas_PDFF_fat_fraction" = tar_var_abb[1],
                              "Liver_PDFF_fat_fraction" = tar_var_abb[2], 
                              "Muscle_fat_infiltration" = tar_var_abb[3],
                              "area_of_pericardial_fat" = tar_var_abb[4])))
cor_results_1[1:2] <- apply(cor_results_1[,1:2],2,function(data){
  ifelse(!data %in% tar_var_abb,toupper(data),data)
})


Item_model <- "Adjust_nor_res_nm" 

for (i in c(1:length(tar_var))){
  
  if(Item_model == "demo_adjust"){
    load(paste0("Process_Data/Fat_Dis/Fat_Causes_Res/",tar_var[i],"_7_100_kocmi.RData"))
    file_save <- "Figures/Fat/Path_Vis/Demo_RM/"
  } else if(Item_model == "no_adjust") {
    load(paste0("Process_Data/Fat_Dis/Fat_Causes/",tar_var[i],"_7_100_kocmi.RData"))
    file_save <- "Figures/Fat/Path_Vis/Demo/"
  } else if(Item_model == "no_adjust_nor") {
    load(paste0("Process_Data/Fat_Dis/Fat_Causes_Nor/",tar_var[i],"_7_100_kocmi.RData"))
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor/"    
  } else if(Item_model == "Adjust_nor_res"){
    load(paste0("Process_Data/Fat_Dis/Fat_Causes_Nor_Res/",tar_var[i],"_7_100_kocmi.RData"))
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor_Res/"   
  } else {
    load(paste0("Process_Data/Fat_Dis/Fat_Causes_Res_NM/",tar_var[i],"_7_100_kocmi_nm.RData"))
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor_Res_NM/"       
  }
  
   res_kocmi <- mutate(res_kocmi,
                      regulator = ifelse(regulator == "target",tar_var_abb[i],toupper(regulator)),
                      target = ifelse(target == "target",tar_var_abb[i],toupper(target)))
   res_kocmi_1 <- merge(res_kocmi,cor_results_1, by = c("regulator","target"))
   res_kocmi_1$t1 <- NULL;res_kocmi_1 <- rename(res_kocmi_1,t1 = correlation)
   
   path_res <- try(
     analyze_pathways(
       weightNet = res_kocmi_1,
       outcome_var = tar_var_abb[i],
       p_threshold = 0.05,cs_threshold = 0.1,
       path_strength_threshold = 0.05), silent = TRUE)
   
   if (inherits(path_res, "try-error")) {
     message(paste("Skip", tar_var_abb[i], "- error"))
     save(res_kocmi_1,file = paste0(file_save,tar_var[i],"_pathways_RData"))
   } else {
     message(paste("Finished", tar_var_abb[i]))
     save(path_res,res_kocmi_1,file = paste0(file_save,tar_var[i],"_pathways_RData"))
  }
} 


cmi_res <- data.frame(); nw_list <- list(); mer_list <- list();

for ( i in c(1:length(tar_var))){

  if(Item_model == "demo_adjust"){
    file_save <- "Figures/Fat/Path_Vis/Demo_RM/"
  } else if(Item_model == "no_adjust"){
    file_save <- "Figures/Fat/Path_Vis/Demo/"
  } else if(Item_model == "no_adjust_nor") {
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor/"
  } else if(Item_model == "no_adjust_nor_res"){
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor_Res/"
  } else {
    file_save <- "Figures/Fat/Path_Vis/Demo_Nor_Res_NM/"
  }
  
  add_new <-gsub("/","",gsub("Figures/Fat/Path_Vis/","",file_save))
  load(file = paste0(file_save,tar_var[i],"_pathways_RData"))
  
  net_stats <- get_network_statistics(weightNet = res_kocmi_1, 
                                      p_adj_threshold = 0.05, 
                                      cs_threshold = 0.1)
  print(net_stats$hub_nodes)
  
  mer_list[[i]] <- main_analysis_vis(weightNet = res_kocmi_1, labels_sub = "",
                                     tar_var_abb[i], 15,show_legend = F)
  
  nw_list[[i]] <- visualize_whole_network(weightNet = res_kocmi_1, 
                                          p_adj_threshold = 0.05, 
                                          cs_threshold = 0.1) +
    labs(title = "", subtitle = "", caption = "") + 
    guides(fill = "none")
  
  cmi_res <- rbind(cmi_res,res_kocmi_1)
  
  vis_plot <- visualize_pathways(analysis_results = path_res)
  pathway_tables <- create_detailed_pathway_tables(path_res)
  saved_results <- save_pathway_analysis(visualization_result = vis_plot,
                                         analysis_results = path_res,
                                         paste0(file_save,tar_var[i]))
  result_summary(weightNet = res_kocmi_1)
  
  save(path_res,res_kocmi_1,file = paste0(file_save,tar_var[i],"_pathways_RData"))
}

network_fin <- plot_grid(plotlist = nw_list,ncol = 1,labels = c(LETTERS[1:4]))
network_fin_u <- plot_grid(plotlist = mer_list,ncol = 2,labels = c(LETTERS[1:4]))

ggsave(filename  =  paste0("Figures/Fat/Direct_Vis/",add_new,"_","full_network.tiff"),
       network_fin,width = 20,height = 35,units = "cm",dpi = 300)
ggsave(filename  =  paste0("Figures/Fat/Direct_Vis/",add_new,"_","Direct_Merge",".tiff"),
       network_fin_u,width = 38,height = 20,units = "cm",dpi = 300)
cmi_res[3:6] <- apply(cmi_res[3:6],2,round,3)
tar_res_cause <- filter(cmi_res,target %in% tar_var_abb & p_adj < 0.05)
write.csv(tar_res_cause,file = paste0("/Users/zhangbing/Documents/Tmp/",add_new,"_cmi_cause.csv"),row.names = F)
tar_res_out <- filter(cmi_res,regulator %in% tar_var_abb & p_adj < 0.05)
write.csv(tar_res_out,file = paste0("/Users/zhangbing/Documents/Tmp/",add_new,"_cmi_outcome.csv"),row.names = F)


filter(cmi_res,target %in% c("Liver PDFF","IGFBP2") & regulator %in% c("Liver PDFF","IGFBP2"))
filter(cor_results_1,target %in% c("Liver PDFF","IGFBP2") & regulator %in% c("Liver PDFF","IGFBP2"))
filter(res_kocmi_1,target %in% c("Liver PDFF","IGFBP2") & regulator %in% c("Liver PDFF","IGFBP2"))
filter(cmi_res,target %in% c("Liver PDFF") & p_adj < 0.05)
filter(cmi_res,target %in% c("Muscle FI") & p_adj < 0.05)
filter(cmi_res,regulator %in% c("Muscle FI") & p_adj < 0.05)


tar_res_u <- filter(cmi_res, p_adj < 0.05 &(regulator %in% tar_var_abb | target %in% tar_var_abb))

save(tar_res_cause,tar_res_out,file = paste0("Process_Data/Fat_Dis/",add_new,"_KOCMI.RData"))




cmi_res_sen <- data.frame(); nw_list_sen <- list();
stat_list_sen <- list(); mer_list_sen <- list()

for ( i in c(1:length(tar_var))){
  load(paste0("Process_Data/Fat_Dis/Fat_Causes_Res/",tar_var[i],"_5100_kocmi.RData"))
  res_kocmi <- mutate(res_kocmi,
                      regulator = ifelse(regulator == "target",tar_var_abb[i],regulator),
                      target = ifelse(target == "target",tar_var_abb[i],target))
  
  detailed_viz <- try(
    visualize_detailed_pathways(
      weightNet = res_kocmi,
      outcome_var = tar_var_abb[i],
      p_threshold = 0.05,
      cs_threshold = 0.1,
      top_direct = 10,
      top_indirect = 15,
      path_strength_threshold = 0.05
    ),
    silent = TRUE
  )
  
  if (inherits(detailed_viz, "try-error")) {
    message(paste("跳过", tar_var_abb[i], "- 函数执行出错"))
  } else {
    message(paste("成功完成", tar_var_abb[i]))
    save(detailed_viz,file = paste0("Figures/Fat/Direct_Indirect/",tar_var[i],"_pathways_RData"))
    print(detailed_viz$combined_plot)
    pathway_tables <- create_detailed_pathway_tables(detailed_viz)
    saved_results <- save_pathway_analysis(detailed_viz, 
                                           paste0("Figures/Fat/Direct_Indirect/",tar_var[i]))
  }

  net_stats <- get_network_statistics(weightNet = res_kocmi, 
                                      p_adj_threshold = 0.05, 
                                      cs_threshold = 0.1)
  print(net_stats$hub_nodes)
  
  comprehensive_plot <- create_network_analysis_report(weightNet = res_kocmi,tar_var_abb[i], 15)
  ggsave(filename  =  paste0("Figures/Fat/Network/Sen/",tar_var[i],".tiff"),comprehensive_plot,
         width = 26,height = 14,units = "cm",dpi = 300)
  
  mer_list_sen[[i]] <- create_network_analysis_report(weightNet = res_kocmi, labels_sub = "",
                                                  tar_var_abb[i], 15,show_legend = F)
  nw_list_sen[[i]] <- plot_causal_network_ggplot(weightNet=res_kocmi, 
                                                 p_adj_threshold = 0.05, 
                                                 cs_threshold = 0.1) +
    labs(title = "", subtitle = "", caption = "") + guides(fill = "none")
  cmi_res_sen <- rbind(cmi_res_sen,res_kocmi)
  ggsave(filename  =  paste0("Figures/Fat/Network/Sen/",tar_var[i],".tiff"), nw_list_sen[[i]],
         width = 22,height = 14,units = "cm",dpi = 300)  
}

network_fin_sen <- plot_grid(plotlist = nw_list_sen,ncol = 2,labels = c(LETTERS[1:4]))

ggsave(filename  =  paste0("Figures/Fat/Network/Sen/",tar_var[i],".tiff"),
       network_fin_sen,width = 26,height = 14,units = "cm",dpi = 300)

network_fin_sen_u <- plot_grid(plotlist = mer_list_sen,ncol = 2,labels = c(LETTERS[1:4]))

ggsave(filename  =  paste0("Figures/Fat/Network/Sen/","merge",".tiff"),
       network_fin_sen_u,width = 38,height = 20,units = "cm",dpi = 300)

tar_res_cause_sen <- filter(cmi_res_sen,target %in% tar_var_abb & p_adj < 0.05)
write.csv(tar_res_cause_sen,file = "/Users/zhangbing/Documents/Tmp/cmi_cause_sen.csv",row.names = F)
tar_res_out_sen <- filter(cmi_res_sen,regulator %in% tar_var_abb & p_adj < 0.05)
write.csv(tar_res_out_sen,file = "/Users/zhangbing/Documents/Tmp/cmi_outcome_sen.csv",row.names = F)

save(tar_res_cause_sen,tar_res_out_sen,
     file = paste0("Process_Data/Fat_Dis/","Causes_Protein_Sen.RData"))


