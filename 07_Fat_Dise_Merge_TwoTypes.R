rm(list = ls())
# Load required packages
library(dplyr)
library(survival)
library(purrr)
library(broom)
library(ggplot2)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)
library(stringi)
library(forcats)
library(ggrepel)
library(patchwork)
library(scales)

source("Codes/00_Basic_Functions.R")
source("Codes_Fat/00_Fat_BasicFunctions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Disease.RData")
load("Process_Data/Covar_Data_Inst2.RData")
dat_death <- read_tsv("ukbrapr_data/death.tsv")
dat_death <- mutate(dat_death,eid = as.numeric(eid))

tar_var <- c("Muscle_fat_infiltration","Liver_PDFF_fat_fraction",
             "Pancreas_PDFF_fat_fraction","area_of_pericardial_fat")

dat.covar <- dat.ana.1
tar_col <- match(tar_var,names(dat.covar))
dat.covar[,tar_col] <- apply(dat.covar[,tar_col],2,scale)

icd_target <- read.csv("Data/ICD_Summ.csv",header = T)

icd_target_list <- lapply(icd_target$Included.ICD.10.Codes,function(data){
  parts <- stri_split_regex(data, ",")[[1]]
  codes <- unlist(lapply(parts, parse_icd_range))
  return(codes)
})
names(icd_target_list) <- icd_target$Grouped.Category


icd_to_group <- data.frame()

for (group_name in names(icd_target_list)) {
  icd_codes <- icd_target_list[[group_name]]
  group_df <- data.frame(
    icd = icd_codes,
    disease_group = group_name,
    stringsAsFactors = FALSE)
  icd_to_group <- bind_rows(icd_to_group, group_df)
}


loop_var <- unique(icd_to_group$disease_group)

dis_reg_list <- lapply(loop_var, function(dis_name){
  
  print(match(dis_name,loop_var))
  dis_name_tem <- filter(icd_to_group,disease_group == dis_name)
  dis_icd_group <- filter(dis_icd,icd %in% as.character(dis_name_tem$icd))
  
  dis_icd_tem <- dis_icd_group %>%
    group_by(eid) %>%
    summarise(dis_date = min(date, na.rm = TRUE))
  
  dat.covar %>%
    left_join(dis_icd_tem,by = "eid") %>%
    left_join(dat_death[,c(2,6)],by = "eid") -> dat.ana
  
  dat.ana.1 <- mutate(dat.ana,
                      age = as.numeric(substr(ins_date2,1,4)) - birth_year,
                      bin = ifelse(!is.na(dis_date),1,0),
                      sex = factor(ifelse(sex == "Male",1,0)),
                      df = dplyr::if_else(is.na(dis_date) & !is.na(date_of_death),
                                          date_of_death,dis_date),
                      df = dplyr::if_else(!is.na(dis_date),dis_date,as.Date("2024-10-30")),
                      t_time = as.numeric(as.Date(df) - as.Date(ins_date2))/365.25)
  dat.ana.2 <- filter(dat.ana.1,t_time > 0)
  
  names(dat.ana.2)[match(c("bin"),names(dat.ana.2))] <- c("status")
  dat.ana.3 <- na_remove(data = dat.ana.2,
                         tar_var = c(tar_var,model_2_var,model_3_var[1:6]))
  
  cardmus_thsd <- dat.ana.3 %>%
    group_by(sex) %>%
    summarise(
      thr_car = quantile(area_of_pericardial_fat, 0.75, na.rm = TRUE),
      thr_muc = quantile(Muscle_fat_infiltration, 0.75, na.rm = TRUE),.groups = 'drop')
  
  # Define gender-specific clinical thresholds
  clinical_thresholds <- list(
    # Male thresholds
    Male = list(
      Pancreas_PDFF = log(10.4),  Liver_PDFF = log(5.0),     
      Muscle = as.numeric(cardmus_thsd[2,3]), 
      area_of_peri = as.numeric(cardmus_thsd[2,2])),
    Female = list(
      Pancreas_PDFF = log(10.4), Liver_PDFF = log(5.0),     
      Muscle = as.numeric(cardmus_thsd[1,3]), 
      area_of_peri = as.numeric(cardmus_thsd[1,2])))
  
  dat_classfied <- create_binary_outcomes(data = dat.ana.3,
                                          thresholds = clinical_thresholds)
  
  binary_var <- c("Pancreas_PDFF_binary","Liver_PDFF_binary",
                   "Muscle_fat_binary","Pericardial_fat_binary")
  
  dat_res_2 <- regression_fit(data = dat_classfied,tar_var = binary_var[1],
                              ori_full = T,rem_var = NULL,
                              add_var = binary_var[-1],ref_level = NULL,
                              exp_var_list = list(model_2_var,model_3_var[1:6]),
                              out_full =TRUE)
  dat_res_2 <- mutate(dat_res_2,id = dis_name);
  return(dat_res_2)
})

dis_reg <- bind_rows(dis_reg_list)

dis_reg <- merge(dis_reg,icd_target,by.x = "id",by.y = "Grouped.Category")
View(dat_tem <- filter(dis_reg,var %in% paste0(binary_var,"1") & model_type == 2 & p_ind < 0.05))
write.csv(dat_tem,file = "/Users/zhangbing/Documents/Tmp/Fat_Dis_Sig_Binary_1.csv",row.names = F)
write.csv(dis_reg,file = "/Users/zhangbing/Documents/Tmp/Fat_Dis_Full_Binary_0.csv",row.names = F)

##################################################################################
###################################  Visualization   #############################
##################################################################################

data_clean <- dat_tem %>%
  mutate(
    # Calculate confidence intervals
    ci_lower = coef - 1.96 * se,
    ci_upper = coef + 1.96 * se,
    # Create significance classification
    significance = case_when(
      p_ind < 0.001 ~ "***",
      p_ind < 0.01 ~ "**", 
      p_ind < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    # Simplify variable names
    var_simple = case_when(
      var == "Muscle_fat_infiltration" ~ "Muscle Fat Infiltration",
      var == "Liver_PDFF_fat_fraction" ~ "Liver Fat Fraction",
      var == "Pancreas_PDFF_fat_fraction" ~ "Pancreas Fat Fraction", 
      var == "area_of_pericardial_fat" ~ "Pericardial Fat Area"
    ),
    # Extract main disease categories
    disease_simple = sapply(strsplit(Description.of.Major.Conditions, ","), function(x) x[1]),
    # Create effect direction
    effect_direction = ifelse(coef > 0, "Positive", "Negative")
  ) %>%
  # Process total frequency commas
  mutate(Total_Freq_Numeric = as.numeric(gsub(",", "", Total.Freq)))

# Set color scheme
var_colors <- c("Muscle Fat Infiltration" = "#E41A1C", 
                "Liver Fat Fraction" = "#377EB8",
                "Pancreas Fat Fraction" = "#4DAF4A", 
                "Pericardial Fat Area" = "#984EA3")

effect_colors <- c("Positive" = "#D55E00", "Negative" = "#0072B2")

# 1. Main Forest Plot: Grouped by Disease System
p1 <- ggplot(data_clean, aes(x = coef, y = fct_reorder(id, coef), color = var_simple)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = 0.3, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_text(aes(label = significance, x = ci_upper + 0.05), 
            position = position_dodge(width = 0.8), size = 3, vjust = 0.5) +
  scale_color_manual(values = var_colors) +
  labs(
    #title = "Association Strength Between Fat Distribution Indicators and Various Diseases",
    #subtitle = "Sorted by disease type, error bars represent 95% confidence intervals",
    x = expression(paste("Regression Coefficient (", beta, ")")),
    y = NULL,
    color = "Fat Indicator",
    caption = "*p<0.05, **p<0.01, ***p<0.001"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0)
  )

print(p1)

# 2. Faceted Association Plot by Disease Chapter
p2 <- ggplot(data_clean, aes(x = coef, y = fct_reorder(id, coef), color = var_simple)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # 使用geom_linerange代替geom_segment，效果更好
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), 
                 position = position_dodge(width = 0.8),
                 alpha = 0.7, size = 1) +
  geom_point(aes(size = -log10(p_ind)), alpha = 0.8, 
             position = position_dodge(width = 0.8)) +
  facet_grid(Chapter ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = var_colors) +
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(1, 4),
    breaks = c(1, 2, 3, 5)
  ) +
  labs(
    x = expression(paste("Regression Coefficient (", beta, ")")),
    y = NULL,
    color = "Fat Indicator"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, size = 9),
    panel.spacing = unit(0.5, "lines")
  )
print(p2)

# 3. Bubble Plot: Disease Frequency vs Effect Size
p3 <- ggplot(data_clean, aes(x = Total_Freq_Numeric, y = abs(coef))) +
  geom_point(aes(color = var_simple, size = -log10(p_ind)), alpha = 0.7) +
  geom_smooth(aes(color = var_simple), method = "lm", se = FALSE, size = 0.5) +
  scale_color_manual(values = var_colors) +
  scale_size_continuous(name = "-log10(p-value)") +
  scale_x_log10(labels = comma) +
  labs(
    #title = "Relationship Between Disease Frequency and Effect Size",
    #subtitle = "Point size indicates significance, color indicates fat indicator type",
    x = "Disease Frequency (log scale)",
    y = expression(paste("Regression Coefficient |", beta, "|")),
    color = "Fat Indicator"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# 4. Association Count Statistics by Fat Indicator
assoc_summary <- data_clean %>%
  group_by(var_simple, effect_direction) %>%
  summarise(count = n(), .groups = "drop")

p4 <- ggplot(assoc_summary, aes(x = var_simple, y = count, fill = effect_direction)) +
  geom_col(position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  scale_fill_manual(values = effect_colors) +
  labs(
    #title = "Number of Associations by Fat Indicator",
    #subtitle = "Categorized by association direction",
    x = NULL,
    y = "Number of Associations",
    fill = "Association Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p4)

# 5. Key Disease Network Plot (Strongest Associations)
top_diseases <- data_clean %>%
  group_by(id) %>%
  slice_max(order_by = abs(coef), n = 1) %>%
  ungroup() %>%
  arrange(desc(abs(coef))) %>%
  head(15)

p5 <- ggplot(top_diseases, aes(x = coef, y = fct_reorder(id, coef))) +
  geom_col(aes(fill = var_simple), alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("β=%.2f\np=%.2e", coef, p_ind), 
                x = ifelse(coef > 0, coef + 0.05, coef - 0.05)), 
            size = 2.8, hjust = ifelse(top_diseases$coef > 0, 0, 1)) +
  scale_fill_manual(values = var_colors) +
  labs(
    #title = "Top 15 Strongest Fat-Disease Associations",
    x =  expression(paste("Regression Coefficient (", beta, ")")),
    y = NULL,
    fill = "Fat Indicator"#,
    #caption = "Shows the fat indicator with the largest effect for each disease"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )

print(p5)

# 6. Association Direction Heatmap (by Disease Chapter)
heatmap_data <- data_clean %>%
  group_by(Chapter, var_simple) %>%
  summarise(
    avg_effect = mean(coef),
    n_associations = n(),
    .groups = "drop"
  ) %>%
  complete(Chapter, var_simple, fill = list(avg_effect = 0, n_associations = 0))

p6 <- ggplot(heatmap_data, aes(x = var_simple, y = Chapter)) +
  geom_tile(aes(fill = avg_effect), color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f\n(n=%d)", avg_effect, n_associations)), 
            size = 2.8, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white", 
    high = "#B2182B",
    midpoint = 0,
    name = "Average Effect Size"
  ) +
  labs(
    #title = "Average Effects of Fat Indicators by Disease Chapter",
    #subtitle = "Cells show average effect size and number of associations",
    x = NULL,
    y = "Disease Chapter"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p6)

# 7. Scatter Plot of p-value Distribution vs Effect Size
p7 <- ggplot(data_clean, aes(x = coef, y = -log10(p_ind))) +
  #geom_point(aes(color = var_simple, size = Total_Freq_Numeric), alpha = 0.7) +
  geom_point(aes(color = var_simple), alpha = 0.7,size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  annotate("text", x = min(data_clean$coef), y = -log10(0.05), 
           label = "p = 0.05", vjust = -0.5, hjust = 0, color = "red", size = 3) +
  scale_color_manual(values = var_colors) +
  scale_size_continuous(
    name = "Disease Frequency",
    trans = "log10",
    breaks = c(10000, 50000, 100000, 300000)
  ) +
  labs(
    # title = "Relationship Between Effect Size and Significance Level",
    #subtitle = "Point size indicates disease frequency, red dashed line indicates p=0.05 threshold",
    x =  expression(paste("Regression Coefficient (", beta, ")")),
    y = "-log10(p-value)",
    color = "Fat Indicator") +
  theme_minimal() + theme(legend.position = c(0.2,0.8)) 

print(p7)


# 8. Combined Graphics Display
combined_plot <- (p2 | p4) / (p5 | p7) +
  plot_annotation(
    title = "Comprehensive Analysis of Fat Distribution and Disease Associations",
    subtitle = "Multi-dimensional visualization display",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )


combined_plot <- plot_grid(p2,plot_grid(p4,p7,ncol = 1),ncol = 2,
                           rel_widths = c(0.7,0.3))

print(combined_plot)

# Save main graphics
ggsave("Figures/Fat/Fat_Dis/fat_disease_main_forest.pdf", p1, width = 14, height = 12, dpi = 300)
ggsave("Figures/Fat/Fat_Dis/fat_disease_heatmap.pdf", p6, width = 10, height = 8, dpi = 300)
ggsave(filename = "Figures/Fat/Fat_Dis/fat_disease_combined.tiff",
       combined_plot, width = 30, height = 18, dpi = 300,units = "cm")

# Generate statistical summary
cat("=== Fat-Disease Association Statistical Summary ===\n")
cat("Total number of associations:", nrow(data_clean), "\n")
cat("Significant associations (p < 0.05):", sum(data_clean$p_ind < 0.05), "\n")
cat("Highly significant associations (p < 0.001):", sum(data_clean$p_ind < 0.001), "\n\n")

cat("Number of associations by fat indicator:\n")
print(table(data_clean$var_simple))

cat("\nAssociation direction distribution:\n")
print(table(data_clean$var_simple, data_clean$effect_direction))

