
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(ggpubr)

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

source("Codes_Fat/00_ML_Shared_functions.R")

sharp_fin <- data.frame(); r2_fin <- data.frame(); pro_name_sel <- data.frame()
sharp_list <- list(); shap_dir_list <- list()

for (i in c(1:length(tar_var))) {
  
   load(file = paste0("Process_Data/Fat_Dis/Continues/", tar_var[i], "_Hierarchical_Models.RData"))

   dat_shap_tem <- all_models_results[[1]]$shap_analysis
   pro_name_dat <- data.frame(pro_name = all_models_results[[1]]$variable_list,
                              var = tar_var[i])
   pro_name_sel <- rbind(pro_name_sel,pro_name_dat)
   sharp_fig_tem <- plot_shap_summary_with_direction(shap_results = dat_shap_tem,
                                                    model_name = "average",
                                                    top_num = 10)
   shap_model <- extract_all_importance(dat_shap_tem)$combined
   shap_dir_list[[i]] <- create_direction_heatmap(shap_model) + ggtitle("") + guides(fill = "none")
  #plot_shap_beeswarm_detailed(dat_shap_tem,model_name = "xgb")
  
  sharp_list[[i]] <- sharp_fig_tem + guides(fill = "none")
   
  sharp_res_list <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_sharp <- dat_tem$shap_analysis
    #plot_shap_summary_with_direction(dat_sharp,model_name = "Model")
    #plot_shap_beeswarm_detailed(dat_sharp,model_name = "xgb")
    dat_sharp <- mutate(dat_sharp$average,model_type = id_row,var = tar_var[i])
    return(dat_sharp)
  })
  sharp_dat <- bind_rows(sharp_res_list);sharp_fin <- rbind(sharp_fin,sharp_dat)

 r2_list <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_r2 <- dat_tem$model_performance; row.names(dat_r2) <- NULL
    dat_r2 <- mutate(dat_r2,model_type = id_row,var = tar_var[i])
    return(dat_r2)
  })
  r2_dat <- bind_rows(r2_list); r2_fin <- rbind(r2_fin,r2_dat)
}  
write.csv(sharp_fin,file = "/Users/zhangbing/Documents/Tmp/sharp.csv",row.names = F)
write.csv(r2_fin,file = "/Users/zhangbing/Documents/Tmp/r2.csv",row.names = F)



sharp_fin_rm <- data.frame(); r2_fin_rm <- data.frame(); pro_name_rm <- data.frame()
sharp_rm_list <- list(); shap_dir_rm_list <- list()

for (i in c(1:length(tar_var))) {
  
  load(file = paste0("Process_Data/Fat_Dis/Continues_Res/", tar_var[i], "_Hierarchical_Models.RData"))
  
  dat_shap_tem <- all_models_results[[1]]$shap_analysis
  
  pro_name_dat <- data.frame(pro_name = all_models_results[[1]]$variable_list,
                             var = tar_var[i])
  pro_name_rm <- rbind(pro_name_rm,pro_name_dat)
  
  sharp_fig_tem <- plot_shap_summary_with_direction(shap_results = dat_shap_tem,
                                                    model_name = "average",top_num = 10)
  shap_model <- extract_all_importance(dat_shap_tem)$combined
  shap_dir_rm_list[[i]] <- create_direction_heatmap(shap_model) + ggtitle("") + guides(fill = "none")
  #plot_shap_beeswarm_detailed(dat_shap_tem,model_name = "xgb")
  
  sharp_rm_list[[i]] <- sharp_fig_tem + guides(fill = "none")
  
  sharp_res_list <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_sharp <- dat_tem$shap_analysis
    #plot_shap_summary_with_direction(dat_sharp,model_name = "Model")
    #plot_shap_beeswarm_detailed(dat_sharp,model_name = "xgb")
    dat_sharp <- mutate(dat_sharp$average,model_type = id_row,var = tar_var[i])
    return(dat_sharp)
  })
  sharp_dat <- bind_rows(sharp_res_list);
  sharp_fin_rm <- rbind(sharp_fin_rm,sharp_dat)
  
  r2_list <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_r2 <- dat_tem$model_performance; row.names(dat_r2) <- NULL
    dat_r2 <- mutate(dat_r2,model_type = id_row,var = tar_var[i])
    return(dat_r2)
  })
  r2_dat <- bind_rows(r2_list); r2_fin_rm <- rbind(r2_fin_rm,r2_dat)
}  
write.csv(r2_fin_rm,file = "/Users/zhangbing/Documents/Tmp/r2_res.csv",row.names = F)
pro_name_sel <- mutate(pro_name_sel,type = "ori");
pro_name_rm <- mutate(pro_name_rm,type = "rm");
pro_name_merge <- rbind(pro_name_sel,pro_name_rm)
write.csv(pro_name_merge,file = "/Users/zhangbing/Documents/Tmp/pro_inter.csv",row.names = F)

data.frame(summarise(group_by(pro_name_merge,var,pro_name),count = n()))


################################################################################
#####################     Coefficient of Determination     #####################
################################################################################

r2_fin_1  <- plot_manu(data = r2_fin)

fig_r2 <- ggplot(r2_fin_1,aes(x = Model, y = R_squared,
                              group = factor(model_type),
                              fill = factor(model_type))) +  
  geom_col(position = "dodge", alpha = 0.9)+
  scale_fill_manual(values = custom_colors[c(3,7,6)]) +
  facet_wrap(~var) + theme_bw(base_size = 12,base_family = "serif") +
  #guides(fill = "none") +
  guides(fill = guide_legend(title = "Model types:")) +
  ylab("Coefficient of Determination") + xlab("Machine learning algorithms") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom")
# 
# ggsave(filename = "Figures/Fat/Prediction/Fig_2_Coefficient_of_Determination.tiff",
#        fig_r2,width = 18,height = 16,units = "cm",dpi = 300)

################################################################################
#####################     Root Mean Square Error     ###########################
################################################################################

fig_rmse <- ggplot(r2_fin_1,aes(x = Model, y = RMSE,
                                group = factor(model_type),
                                fill = factor(model_type))) +  
  geom_col(position = "dodge", alpha = 0.9)+
  scale_fill_manual(values = custom_colors[c(3,7,6)]) +
  facet_wrap(~var,scales = "free_y") +
  theme_bw(base_size = 12,base_family = "serif") +
  guides(fill = guide_legend(title = "Model types")) +
  ylab("Root Mean Square Error") + xlab("Machine learning algorithms") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom")

ggsave(filename = "Figures/Fat/Prediction/Fig_S3_Root_Mean_Square_Error.tiff",
       fig_rmse,width = 24,height = 16,units = "cm",dpi = 300)

################################################################################
#####################    classification visualization   ########################
################################################################################

sharp_fin_clu <- data.frame(); auc_fin <- data.frame(); shap_clu_list <- list()
shap_clu_dir_list <- list();pro_name_clu <- data.frame()

for (i in c(1:length(tar_var))) {
  
  load(file = paste0("Process_Data/Fat_Dis/Classification/", tar_var[i], "_Hierarchical_Classification_Models.RData"))

  ###################### Model performance  ##########################
  
  dat_shap_tem <- all_models_results[[1]]$shap_analysis
  shap_fig_tem <- plot_shap_summary_with_direction(dat_shap_tem,model_name = "average",top_num = 10)
  #plot_shap_beeswarm_detailed(dat_shap_tem$shap_analysis,model_name = "xgb")
  shap_clu_list[[i]] <- shap_fig_tem + guides(fill = "none")

  pro_name_dat <- data.frame(pro_name = all_models_results[[1]]$variable_list,
                             var = tar_var[i])
  pro_name_clu <- rbind(pro_name_clu,pro_name_dat)
  
  shap_clu_model <- extract_all_importance(dat_shap_tem)$combined
  shap_clu_dir_list[[i]] <- create_direction_heatmap(shap_clu_model) + ggtitle("") + guides(fill = "none")
  
  sharp_list_tem <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_sharp <- dat_tem$shap_analysis
    dat_sharp <- mutate(dat_sharp$average,model_type = id_row,var = tar_var[i])
    return(dat_sharp)
  })
  sharp_dat <- bind_rows(sharp_list_tem);sharp_fin_clu <- rbind(sharp_fin_clu,sharp_dat)

  ###################### AUC  ##########################
  
  auc_list <- lapply(1:length(all_models_results),function(id_row){
    dat_tem <- all_models_results[[id_row]]
    dat_auc <- dat_tem$model_performance; row.names(dat_auc) <- NULL
    dat_auc <- mutate(dat_auc,model_type = id_row,var = tar_var[i])
    return(dat_auc)
  })
  auc_dat <- bind_rows(auc_list); auc_fin <- rbind(auc_fin,auc_dat)
}  
write.csv(auc_fin,file = "/Users/zhangbing/Documents/Tmp/auc.csv",row.names = F)

sharp_fin %>%
  filter(model_type == 1) %>%
  group_by(var) %>%
  arrange(desc(Mean_Abs_SHAP)) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything()) %>%
  ungroup()  %>% 
  mutate(type = "continuous")  -> sharp_fin_10

sharp_fin_clu %>%
  filter(model_type == 1) %>%
  group_by(var) %>%
  arrange(desc(Mean_Abs_SHAP)) %>%
  slice_head(n = 10) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything()) %>%
  ungroup()  %>% 
  mutate(type = "binary") -> sharp_fin_clu_10

sharp_fin_fin <- rbind(sharp_fin_clu_10,sharp_fin_10)
  
write.csv(sharp_fin_fin[,-c(5:7,12)],file = "/Users/zhangbing/Documents/Tmp/shap_fin.csv",row.names = F)


################################################################################
#####################     importance visualization         #####################
################################################################################

shap_fig <- plot_grid(plotlist = sharp_list,ncol = 2,labels = LETTERS[1:4])
ggsave(filename = "Figures/Fat/Prediction/Fig_S4_Improtance.tiff",shap_fig,
       width = 24,height = 18,units = "cm",dpi = 300)
shap_fig_u <- plot_grid(plotlist = sharp_list,ncol = 1,labels = LETTERS[1:4])

shap_clu_fig <- plot_grid(plotlist = shap_clu_list,ncol = 2,labels = LETTERS[1:4])
shap_clu_fig_u <- plot_grid(plotlist = shap_clu_list,ncol = 1,labels = "")
ggsave(filename = "Figures/Fat/Prediction/Fig_S5_Improtance.tiff",shap_fig,
       width = 24,height = 18,units = "cm",dpi = 300)

fig_sharp_fin <- plot_grid(shap_fig_u,shap_clu_fig_u,ncol = 2)

ggsave(filename = "Figures/Fat/Prediction/Fig_3_Improtance.tiff",fig_sharp_fin,
       width = 18,height = 20,units = "cm",dpi = 300)


################################################################################
#####################     AUC and kappa visualization         ##################
################################################################################

auc_fin_1 <- plot_manu(auc_fin)

fig_AUC <- ggplot(auc_fin_1,aes(x = Model, y = AUC,
                                group = factor(model_type),
                                fill = factor(model_type))) +  
  geom_col(position = "dodge", alpha = 0.9)+
  scale_fill_manual(values = custom_colors[c(3,7,6)]) +
  facet_wrap(~var) +
  scale_y_continuous(limits = c(0,1))+
  theme_bw(base_size = 12,base_family = "serif") +
  guides(fill = guide_legend(title = "Model types")) +
  ylab("Area Under the ROC Curve") + xlab("Machine learning algorithms") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom")
  
fig_kappa <- ggplot(auc_fin_1,aes(x = Model, y = Kappa,
                                group = factor(model_type),
                                fill = factor(model_type))) +  
  geom_col(position = "dodge", alpha = 0.9)+
  scale_fill_manual(values = custom_colors[c(3,7,6)]) +
  facet_wrap(~var) +
  scale_y_continuous(limits = c(0,1))+
  theme_bw(base_size = 12,base_family = "serif") +
  guides(fill = guide_legend(title = "Model types")) +
  ylab("Cohen’s Kappa Coefficient") + xlab("Machine learning algorithms") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom")

ggsave(filename = "Figures/Fat/Prediction/Fig_S4_Improtance_Kappa.tiff",fig_kappa,
       width = 26,height = 18,units = "cm",dpi = 300)

fig_merge <- plot_grid(fig_r2+guides(fill = "none"),fig_AUC,ncol = 1,
                       labels = c("A","B"))

ggsave(filename = "Figures/Fat/Prediction/Fig_2_Predication_Score.tiff",fig_merge,
       width = 28,height = 20,units = "cm",dpi = 300)


################################################################################
#####################     direction visualization          #####################
################################################################################

rel_width_list <- list(c(0.5,0.5),c(0.65,0.35),c(0.65,0.35),c(0.35,0.65))

for (i in c(1:length(tar_var))){
  
  fig_tem <- plot_grid(shap_dir_list[[i]] + theme(
    axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()),
    shap_clu_dir_list[[i]],ncol = 1,rel_heights = rel_width_list[[i]],
    labels = c("A","B"))
  ggsave(filename  =  paste0("Figures/Fat/Shap_Direct/",tar_var[i],".tiff"),fig_tem,
         width = 22,height = 24,units = "cm",dpi = 300)
}



pro_name_sel <- mutate(pro_name_sel,type = "ori");


fig_pro_list <- list(); model_res <- data.frame()

for (i in tar_var) {
  
  continuous_proteins <- toupper(filter(pro_name_sel, var == i)[, 1])
  binary_proteins <- toupper(filter(pro_name_clu, var == i)[, 1])
  all_proteins <- unique(c(continuous_proteins, binary_proteins))
  
  # Create protein classification dataframe
  protein_categories <- data.frame(
    Protein = all_proteins,
    Category = case_when(
      all_proteins %in% intersect(continuous_proteins, binary_proteins) ~ "Both",
      all_proteins %in% continuous_proteins & !(all_proteins %in% binary_proteins) ~ "Continuous Only",
      all_proteins %in% binary_proteins & !(all_proteins %in% continuous_proteins) ~ "Binary Only",
      TRUE ~ "Other"
    ), type = i)
  model_res <- rbind(model_res,protein_categories)
  
  # Prepare data for network plot
  adj_matrix <- matrix(0, nrow = length(all_proteins), ncol = 2)
  rownames(adj_matrix) <- all_proteins
  colnames(adj_matrix) <- c("Continuous", "Binary")
  
  adj_matrix[all_proteins %in% continuous_proteins, "Continuous"] <- 1
  adj_matrix[all_proteins %in% binary_proteins, "Binary"] <- 1
  
  network_data <- as.data.frame(adj_matrix) %>%
    mutate(Protein = rownames(.)) %>%
    pivot_longer(cols = c(Continuous, Binary), 
                 names_to = "Model_Type", 
                 values_to = "Selected") %>%
    filter(Selected == 1) %>%
    left_join(protein_categories, by = "Protein") 
  
  # Create color encoding: whether in intersection
  network_data <- network_data %>%
    mutate(
      In_Intersection = ifelse(Category == "Both", "In Intersection", "Not In Intersection"),
      Model_Type = factor(Model_Type,levels = c("Continuous","Binary"))
    )
  
  # Create network plot (without legend)
  network_plot <- ggplot(network_data, aes(x = Model_Type, y = Protein)) +
    geom_line(aes(group = Protein, color = In_Intersection), 
              alpha = 0.8, linewidth = 1.2) +
    geom_point(aes(color = In_Intersection, shape = Category), size = 2.5) +
    scale_color_manual(
      name = "Intersection Status",
      values = c(
        "In Intersection" = pal_npg()(2)[1],
        "Not In Intersection" = pal_npg()(2)[2]
      ),
      labels = c("In Intersection" = "Selected in Both Models", 
                 "Not In Intersection" = "Selected in One Model")
    ) +
    scale_shape_manual(
      name = "Selection Pattern",
      values = c(
        "Both" = 18,
        "Continuous Only" = 16,
        "Binary Only" = 17),
      labels = c(
        "Both" = "Both Models",
        "Continuous Only" = "Continuous Only",
        "Binary Only" = "Binary Only"
      )) +
    labs(
      x = "Model Type",
      y = "Protein") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 9),
      panel.grid.major.x = element_blank(),
      legend.position = "none",  # Remove individual plot legend
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
    )
  
  fig_pro_list[[match(i, tar_var)]] <- network_plot
}

group_dat <- data.frame(summarise(group_by(model_res,type,Category),count = n()))
all_dat <- data.frame(summarise(group_by(model_res,type),count = n()))
group_dat_1 <- merge(group_dat,all_dat,by = "type")
group_dat_1 <- mutate(group_dat_1,ratio = count.x/count.y)


get_legend <- function(my_ggplot) {
  tmp <- ggplot_gtable(ggplot_build(my_ggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Create a plot with legend to extract it
legend_plot <- ggplot(network_data, aes(x = Model_Type, y = Protein)) +
  geom_line(aes(group = Protein, color = In_Intersection), 
            alpha = 0.8, linewidth = 1.2) +
  geom_point(aes(color = In_Intersection, shape = Category), size = 2.5) +
  scale_color_manual(
    name = "Intersection Status",
    values = c(
      "In Intersection" = pal_npg()(2)[1],
      "Not In Intersection" = pal_npg()(2)[2]
    ),
    labels = c("In Intersection" = "Selected in Both Models", 
               "Not In Intersection" = "Selected in One Model")) +
  scale_shape_manual(
    name = "Selection Pattern",
    values = c(
      "Both" = 18,
      "Continuous Only" = 16,
      "Binary Only" = 17),
    labels = c(
      "Both" = "Both Models",
      "Continuous Only" = "Continuous Only",
      "Binary Only" = "Binary Only"
    )) +
  theme(legend.position = "bottom", legend.box = "vertical")

shared_legend <- get_legend(legend_plot)

# Method 1: Merge four plots into a grid using plot_grid
combined_plot <- plot_grid(
  plotlist = fig_pro_list,
  nrow = 1,
  ncol = 4,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  align = "hv",
  axis = "lrtb")

# Add shared legend to the bottom of the combined plot
final_plot <- plot_grid(
  combined_plot,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.1)  # Adjust height ratio of legend
)

# Display final plot
print(final_plot)

ggsave("Figures/Fat/Prediction/Fig_S3_lasso_proteins_network.png", 
       plot = final_plot,width = 12, height = 12, dpi = 300)


pro_name <- read.csv("Data/Protein/UKB_PPP.csv",header = T)
pro_name <- mutate(pro_name,meaning = tolower(meaning))
  
var_imp_u <- merge(var_imp,pro_name[,-1],by.x = "Protein",by.y = "meaning",all.x = T)
var_imp_u <- mutate(var_imp_u,dir_ind = ifelse(coef < 0,0,1))


  
  