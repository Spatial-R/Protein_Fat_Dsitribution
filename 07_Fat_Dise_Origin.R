
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(survival)
library(readr)
library(survminer)

source("Codes/00_Basic_Functions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Disease.RData")
load("Process_Data/Covar_Data_Inst2.RData")
dat_death <- read_tsv("ukbrapr_data/death.tsv")
dat_death <- mutate(dat_death,eid = as.numeric(eid))

dat.covar <- dat.ana.1

dis_sum <- data.frame(table(dis_icd$icd))
#write.csv(dis_sum,file = "/Users/zhangbing/Documents/Tmp/ICD10.csv",row.names = F)

dis_sum <- filter(dis_sum,Freq > 1000)

dat_body_fin <- filter(dat_body_fin_u, instance == 2)

tar_var <- c("Muscle_fat_infiltration","Liver_PDFF_fat_fraction",
             "Pancreas_PDFF_fat_fraction","area_of_pericardial_fat")
dat.covar[,match(tar_var,names(dat.covar))] <- apply(dat.covar[,match(tar_var,names(dat.covar))],
                                                     2,scale)

res_fin <- data.frame()

for (i in c(1:nrow(dis_sum))) {
  
  print(i)
  dat_dem <- filter(dis_icd,icd == dis_sum[i,1])
  names(dat_dem)[2] <- "dis_date"
  
  dat.covar %>%
    left_join(dat_dem[,1:2],by = "eid") %>%
    left_join(dat_death[,c(2,6)],by = "eid") -> dat.ana
    #left_join(dat_cor,by = "eid") -> dat.ana
  
  dat.ana.1 <- mutate(dat.ana,
                      age = as.numeric(substr(ins_date2,1,4)) - birth_year,
                      bin = ifelse(!is.na(dis_date),1,0),
                      df = dplyr::if_else(is.na(dis_date) & !is.na(date_of_death),
                                          date_of_death,dis_date),
                      df = dplyr::if_else(!is.na(dis_date),dis_date,as.Date("2024-10-30")),
                      t_time = as.numeric(as.Date(df) - as.Date(ins_date2))/365.25)
  dat.ana.2 <- filter(dat.ana.1,t_time > 0)
  
  names(dat.ana.2)[match(c("bin"),names(dat.ana.2))] <- c("status")
  dat.ana.3 <- na_remove(data = dat.ana.2,
                         tar_var = c(tar_var,model_2_var,model_3_var[1:6]))
  
  dat_res_2 <- regression_fit(data = dat.ana.3,tar_var = tar_var[1],
                              ori_full = T,rem_var = NULL,
                              add_var = tar_var[-1],ref_level = NULL,
                              exp_var_list = list(model_2_var,model_3_var[1:6]),out_full =TRUE)
  dat_res_2 <- mutate(dat_res_2,id = dis_sum[i,1]);
  res_fin <- rbind(res_fin,dat_res_2)
}

save(res_fin,file = "Process_Data/Fat_Dis/Fat_Diseases.RData")

res_fin_u1 <- filter(res_fin,var %in% tar_var & model_type == 2)
res_fin_u1 <- mutate(res_fin_u1,p_ind_u = ifelse(p_ind < 0.05,1,0))
res_fin_u2 <- merge(res_fin_u1,dis_name_1,by.x = 'id',by.y = "des_1")
res_fin_u3 <- filter(res_fin_u2, p_ind_u == 1) 
write.csv(res_fin_u3[,c(1,6,8,2,4)],file = "Process_Data/Fat_Dis/Fat_Dis_Sig.csv",row.names = F)
write.csv(res_fin_u3[,c(1,6,8,2,4)],file = "/Users/zhangbing/Documents/Tmp/Fat_Dis_Sig.csv",row.names = F)
write.csv(res_fin_u2[,c(6,8,2,4)],file = "Process_Data/Fat_Dis/Fat_Dis_All.csv",row.names = F)


# 加载必要的包
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(patchwork)
library(ggrepel)
library(gridExtra)
library(forcats)

# Set theme for scientific publication
theme_sci <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey80"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.2, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "lines")
    )
}

# Set global theme
theme_set(theme_sci())

# Read and preprocess data
#fat_data <- read.csv("Fat_Dis_Sig.csv", stringsAsFactors = FALSE)



# Set scientific theme
theme_sci <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey80"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")
    )
}

theme_set(theme_sci())

# Color palette
fat_colors <- c(
  "Muscle Fat\nInfiltration" = "#1f77b4",
  "Liver PDFF\nFat Fraction" = "#ff7f0e", 
  "Pancreas PDFF\nFat Fraction" = "#2ca02c",
  "Pericardial\nFat Area" = "#d62728"
)

# Function to create detailed visualization for a specific disease category
create_disease_category_plot <- function(disease_cat, data) {
  
  # Filter data for the specific disease category
  category_data <- data %>%
    filter(disease_category == disease_cat) %>%
    mutate(
      disease_short = str_trunc(dis, width = 50),
      significance = case_when(
        p_ind < 0.001 ~ "***",
        p_ind < 0.01 ~ "**", 
        p_ind < 0.05 ~ "*",
        TRUE ~ "NS"
      )
    )
  
  if (nrow(category_data) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available") +
             labs(title = disease_cat) +
             theme_void())
  }
  
  # 1. Association coefficient plot
  p1 <- ggplot(category_data, 
               aes(x = fct_reorder(disease_short, coef), 
                   y = coef, fill = fat_variable)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
    geom_text(aes(label = significance, y = ifelse(coef >= 0, coef + 0.05, coef - 0.05)),
              position = position_dodge(width = 0.7), size = 3) +
    coord_flip() +
    scale_fill_manual(values = fat_colors, name = "Fat Variable") +
    labs(
      title = paste("Association Coefficients -", disease_cat),
      subtitle = "Significance: *P<0.05, **P<0.01, ***P<0.001",
      x = "Disease",
      y = "Association Coefficient (β)"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    )
  
  # 2. Volcano plot for the category
  p2 <- ggplot(category_data, aes(x = coef, y = neg_log10_p, color = fat_variable)) +
    geom_point(aes(size = abs(coef)), alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(aes(label = str_trunc(dis, 30)), 
                    size = 2.5, max.overlaps = 10, segment.size = 0.2) +
    scale_color_manual(values = fat_colors, name = "Fat Variable") +
    scale_size_continuous(range = c(2, 6), name = "|β|") +
    labs(
      title = paste("Volcano Plot -", disease_cat),
      x = "Association Coefficient (β)",
      y = "-log₁₀(P-value)"
    )
  
  # 3. Fat variable distribution within category
  p3 <- category_data %>%
    count(fat_variable, significance) %>%
    ggplot(aes(x = fat_variable, y = n, fill = significance)) +
    geom_col(position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("***" = "#003366", "**" = "#336699", "*" = "#6699cc", "NS" = "#cccccc"),
                      name = "Significance") +
    labs(
      #title = paste("Association Count by Fat Variable -", disease_cat),
      x = "Fat Variable",
      y = "Number of Associations"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 4. Summary statistics table (as plot) - using ggpubr
  # 首先，创建摘要统计表
  summary_stats <- category_data %>%
    group_by(fat_variable) %>%
    summarise(
      n = n(),
      n_sig = sum(p_ind < 0.05),
      mean_beta = round(mean(coef), 3),
      sd_beta = round(sd(coef), 3),
      min_beta = round(min(coef), 3),
      max_beta = round(max(coef), 3),
      .groups = 'drop'
    )
  
  # 将摘要统计表转换为适合ggtexttable的格式
  table_data <- summary_stats %>%
    rename(
      `Fat Variable` = fat_variable,
      `N` = n,
      `Sig` = n_sig,
      `Mean β` = mean_beta,
      `SD β` = sd_beta,
      `Min β` = min_beta,
      `Max β` = max_beta
    )
  
  # 使用ggpubr创建表格图
  p4 <- ggpubr::ggtexttable(table_data, rows = NULL, 
                            theme = ggpubr::ttheme("classic")) +
    labs(title = paste("Summary Statistics -", disease_cat)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = paste('Detailed Analysis:', disease_cat),
                    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))
  
  return(combined_plot)
}

# Read and preprocess data (using the same preprocessing as before)
fat_data <- res_fin_u3[,c(1,6,8,2,4)]

fat_data_clean <- fat_data %>%
  mutate(
    fat_variable = case_when(
      var == "Muscle_fat_infiltration" ~ "Muscle Fat\nInfiltration",
      var == "Liver_PDFF_fat_fraction" ~ "Liver PDFF\nFat Fraction", 
      var == "Pancreas_PDFF_fat_fraction" ~ "Pancreas PDFF\nFat Fraction",
      var == "area_of_pericardial_fat" ~ "Pericardial\nFat Area",
      TRUE ~ var
    ),
    disease_category = case_when(
      str_detect(id, "^A|^B") ~ "Infectious Diseases",
      str_detect(id, "^D") ~ "Hematologic Diseases",
      str_detect(id, "^E") ~ "Metabolic Diseases",
      str_detect(id, "^F") ~ "Psychiatric Disorders",
      str_detect(id, "^G") ~ "Neurological Diseases",
      str_detect(id, "^H") ~ "Ophthalmic Diseases",
      str_detect(id, "^I") ~ "Cardiovascular Diseases",
      str_detect(id, "^J") ~ "Respiratory Diseases",
      str_detect(id, "^K") ~ "Digestive Diseases",
      str_detect(id, "^L") ~ "Dermatologic Diseases",
      str_detect(id, "^M") ~ "Musculoskeletal Diseases",
      str_detect(id, "^N") ~ "Genitourinary Diseases",
      TRUE ~ "Other Diseases"
    ),
    neg_log10_p = -log10(p_ind)
  )

# Create detailed visualizations for each category
cat1_plot <- create_disease_category_plot("Cardiovascular Diseases", fat_data_clean)
cat2_plot <- create_disease_category_plot("Digestive Diseases", fat_data_clean) 
cat3_plot <- create_disease_category_plot("Metabolic Diseases", fat_data_clean)
cat4_plot <- create_disease_category_plot("Musculoskeletal Diseases", fat_data_clean)
cat5_plot <- create_disease_category_plot("Neurological Diseases", fat_data_clean)

# Display plots
print(cat1_plot)
print(cat2_plot)
print(cat3_plot)
print(cat4_plot)

# Save plots
ggsave("Figures/Fat/Fat_Dis/Cardiovascular_Diseases_Detailed.tiff", cat1_plot, 
       width = 25, height = 20, units = "cm", dpi = 300, compression = "lzw")
ggsave("Figures/Fat/Fat_Dis/Digestive_Diseases_Detailed.tiff", cat2_plot, 
       width = 25, height = 20, units = "cm", dpi = 300, compression = "lzw")
ggsave("Figures/Fat/Fat_Dis/Metabolic_Diseases_Detailed.tiff", cat3_plot, 
       width = 25, height = 20, units = "cm", dpi = 300, compression = "lzw")
ggsave("Figures/Fat/Fat_Dis/Musculoskeletal_Diseases_Detailed.tiff", cat4_plot, 
       width = 25, height = 20, units = "cm", dpi = 300, compression = "lzw")

# Additional: Comparative analysis across the four categories
comparison_data <- fat_data_clean %>%
  filter(disease_category %in% c("Cardiovascular Diseases", "Digestive Diseases", 
                                 "Metabolic Diseases", "Musculoskeletal Diseases"))

# Comparative summary statistics
comparison_summary <- comparison_data %>%
  group_by(disease_category, fat_variable) %>%
  summarise(
    n_associations = n(),
    n_significant = sum(p_ind < 0.05),
    percent_significant = round(n_significant / n_associations * 100, 1),
    mean_coef = round(mean(coef), 3),
    median_coef = round(median(coef), 3),
    .groups = 'drop'
  )

# Comparative heatmap
comp_heatmap <- ggplot(comparison_summary, 
                       aes(x = fat_variable, y = disease_category, fill = mean_coef)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f\n(%d/%d)", mean_coef, n_significant, n_associations)), 
            size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "Mean β") +
  labs(
    #title = "Comparative Analysis: Mean Association Coefficients",
    subtitle = "Format: Mean β (Significant/Total Associations)",
    x = "Fat Distribution Variable",
    y = "Disease Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Comparative bar plot: percentage of significant associations
comp_bar <- comparison_summary %>%
  ggplot(aes(x = disease_category, y = percent_significant, fill = fat_variable)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = fat_colors, name = "Fat Variable") +
  labs(
    #title = "Percentage of Statistically Significant Associations",
    x = "Disease Category",
    y = "Percentage Significant (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine comparative plots
comparative_plot <- comp_heatmap / comp_bar +
  plot_annotation(title = "Cross-Category Comparative Analysis",
                  theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))

print(comparative_plot)
ggsave("Four_Categories_Comparative.tiff", comparative_plot, 
       width = 20, height = 16, units = "cm", dpi = 300, compression = "lzw")

# Print summary statistics for manuscript
cat("\n=== COMPARATIVE SUMMARY STATISTICS ===\n\n")
print(comparison_summary)

# Calculate overall statistics for each category
category_overview <- comparison_data %>%
  group_by(disease_category) %>%
  summarise(
    total_associations = n(),
    total_significant = sum(p_ind < 0.05),
    percent_significant = round(total_significant / total_associations * 100, 1),
    mean_absolute_coef = round(mean(abs(coef)), 3),
    strongest_positive = round(max(coef), 3),
    strongest_negative = round(min(coef), 3),
    .groups = 'drop'
  )

cat("\n=== CATEGORY OVERVIEW ===\n")
print(category_overview)
