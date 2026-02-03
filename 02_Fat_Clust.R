
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(factoextra)
library(ggcorrplot)
library(mclust)
source("Codes_Fat/00_Fat_BasicFunctions.R")

tar_var <- c("Pancreas_PDFF_fat_fraction","Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration","area_of_pericardial_fat")

source("Codes/00_Basic_Functions.R")
source("Codes_Fat/00_Fat_BasicFunctions.R")

load(file = "Process_Data/Fat_Dis/Model_Data.RData")

dat.cluster <- dat.ana.u.1[,match(c(tar_var,"age","sex"),names(dat.ana.u.1))]
dat.cluster <- batch_remove_effects(dat.cluster,tar_var,covariate_vars = c("age","sex") )
feature_sample_tem_0 <- data.frame(scale(dat.cluster[,c(1:4)]))
names(feature_sample_tem_0) <- gsub("_residual","",names(feature_sample_tem_0))

outlier_scores <- apply(feature_sample_tem_0, 2, function(x) {
  abs(x - median(x)) / mad(x)
})
outlier_flag <- apply(outlier_scores, 1, function(x) any(x > 5))
cat("We found ", sum(outlier_flag), "outiers \n")

feature_sample_tem <- feature_sample_tem_0[!outlier_flag, ]
#plot(fviz_nbclust(feature_sample_tem_0, kmeans, method = "silhouette", k.max = 10))
#feature_sample_tem <- residuals_df

mc_res_tem <- Mclust(feature_sample_tem);summary(mc_res_tem)
plot(mc_res_tem, what = "BIC", 
     ylab = "BIC value",xlab = "Number of components",
     main = "Model Selection Based on BIC")
print(mc_res_tem$G)
cluster_num <- mc_res_tem$G
#if(cluster_num > 4) cluster_num <- 4

mc_res_tem_1 <- Mclust(feature_sample_tem,G = cluster_num);

dat_mclust_tem <- data.frame(feature_sample_tem) %>%
  mutate(cluster = as.factor(mc_res_tem_1$classification))

(cluster_props <- table(dat_mclust_tem$cluster)/nrow(dat_mclust_tem))

used_colors <- get_palette("lancet", cluster_num)#[c(2,1,4,3)]

fig_1_tem <- fviz_mclust(mc_res_tem_1, what = "classification",
                         geom = "point", palette = used_colors) +
  ggtitle("")

cluster_means_tem <- dat_mclust_tem[,-6] %>%
  group_by(cluster) %>%
  dplyr::select(where(is.numeric)) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-cluster, names_to = "Variable", values_to = "Mean_Value")

cluster_means_tem <- mutate(cluster_means_tem,
                            Variable = factor(Variable,levels = tar_var,
                                              labels = fat_abbr))

ordered_legend_labels <- paste0("Cluster ", levels(dat_mclust_tem$cluster), 
                                " (", round(cluster_props*100, 1), "%)")


fig_2_tem <- ggplot(cluster_means_tem, aes(x = Variable, y = Mean_Value, 
                                           color = cluster, group = cluster)) +
  geom_line() +
  geom_point(aes(shape = cluster,
                 color = cluster),size = 2) +
  # coord_flip()+
  scale_color_manual(
    name = "Clusters (%)",
    values = used_colors,
    labels = ordered_legend_labels
  ) +
  scale_shape_manual(
    name = "Clusters (%)",  # 使用相同的图例标题
    #values = c(16,17,15,3),
    labels = ordered_legend_labels  # 使用相同的标签
  ) +
  labs(#title = "Cluster Profile Plot",
    y = "Mean Values",
    x = "") +
  theme_classic(base_family = "serif",base_size = 12) +
  theme(#axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.8,0.85))


ggplot(cluster_means_tem, aes(x = Variable, y = Mean_Value, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_color_manual(values = used_colors)+
  coord_flip() + # 翻转坐标轴，让变量名更易读
  labs(title = "Cluster Means by Variable",
       y = "Mean Values",
       x = "") + theme_bw()

fig_fin <- plot_grid(fig_1_tem+guides(shape="none",color = "none",fill = "none"),
                     fig_2_tem,ncol = 2,
                                labels = paste(LETTERS[1:2], sep = ""))

save(fig_1_tem,fig_2_tem,dat_mclust_tem,file = "Figures/Fat/Clust.RData")

ggsave(filename ="Figures/Fat/Fat_Cluster.tiff",fig_fin,
       width = 22,height = 13,units = "cm",dpi = 300)

