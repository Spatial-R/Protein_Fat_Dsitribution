
rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(scales)
library(patchwork)
library(grid)
library(gridExtra)
library(viridis)
library(ggsci)

source("Codes_Fat/05_Protein_Enrich_Functions.R")

# Main analysis workflow
tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

dir_target <- "Figures/Fat/Enrich/"


go_list <- list(); kegg_list <- list(); res_dat <- data.frame()

for (i in 1:4) {
  
  cat("Analyzing:", tar_var[i], "\n")
  
  # Load data
  load(paste0("Figures/Fat/Path_Vis/Demo_Nor_Res_NM/", tar_var[i], "_pathways_RData"))
  
  # Get target proteins
  target_proteins <- unique(c(path_res$direct_edges$regulator, 
                              path_res$indirect_paths$source,
                              unlist(path_res$indirect_paths$mediators)))
  
  res_fin_tem <- run_complete_enrichment_analysis(
                                   protein_names = target_proteins,
                                   output_dir = paste0(dir_target,tar_var[i],"/"),
                                   p_cutoff = 0.05,
                                   q_cutoff = 0.2,
                                   create_plots = TRUE)
  go_list[[i]] <- res_fin_tem$figures$go_figure;
  kegg_list[[i]] <- res_fin_tem$figures$kegg_plot;
  
  combined_enr <- bind_rows(
    BP = data.frame(res_fin_tem$go_results$BP),
    MF = data.frame(res_fin_tem$go_results$MF),  # 或 MP
    CC = data.frame(res_fin_tem$go_results$CC),
    kegg = data.frame(res_fin_tem$kegg_result),
    .id = "Category"
  )
  combined_enr <- mutate(combined_enr,fat = tar_var[i])
  res_dat <- rbind(res_dat,combined_enr)
  save(res_fin_tem, target_proteins,file = paste0(dir_target, tar_var[i],"/", "Results.RData"))
}

write.csv(res_dat,file = "/Users/zhangbing/Documents/Tmp/Enrich.csv",row.names = F)
fig_1 <- plot_grid(plotlist = go_list,ncol = 1,labels = c("A","B","C","D"))
ggsave("Figures/Fat/Enrich/Fig_4.png", 
       plot = fig_1,width = 10, height = 12, dpi = 300)

fig_2 <- plot_grid(kegg_list[[2]],kegg_list[[3]],kegg_list[[4]],ncol = 1,labels = c("A","B","C"))
ggsave("Figures/Fat/Enrich/Fig_S12.png", 
       plot = fig_2,width = 10, height = 7, dpi = 300)

