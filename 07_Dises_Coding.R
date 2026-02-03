

# 创建层次化合并函数，包含病例数阈值
create_hierarchical_groups_with_threshold <- function(icd_data, detailed_groups, threshold = 300) {
  # 创建数据框存储结果
  result <- data.frame()
  
  # 处理每个详细分组
  for (group_name in names(detailed_groups)) {
    codes <- detailed_groups[[group_name]]
    # 从原始数据中提取这些编码
    group_data <- icd_data %>% 
      filter(icd %in% codes) %>%
      summarise(
        Disease_Group = group_name,
        ICD_Codes = paste(codes, collapse = ", "),
        Total_Cases = sum(Freq),
        Number_of_Codes = length(codes),
        Mean_Cases_per_Code = round(mean(Freq), 1)
      )
    result <- bind_rows(result, group_data)
  }
  
  # 计算每个组在总人群中的患病率
  total_population <- sum(icd_data$Freq)
  result <- result %>%
    mutate(
      Prevalence_per_100k = round((Total_Cases / total_population) * 100000, 1),
      Proportion_of_Total = round((Total_Cases / total_population) * 100, 2),
      # 标记是否低于阈值
      Below_Threshold = Total_Cases < threshold
    ) %>%
    arrange(desc(Total_Cases))
  
  return(result)
}

# 应用详细合并
detailed_merged_data <- create_hierarchical_groups_with_threshold(icd_data, detailed_disease_groups, threshold = 300)

# 识别需要合并的低频疾病组
low_frequency_groups <- detailed_merged_data %>%
  filter(Below_Threshold) %>%
  pull(Disease_Group)

cat("需要合并的低频疾病组数量:", length(low_frequency_groups), "\n")
cat("低频疾病组:", paste(low_frequency_groups, collapse = ", "), "\n")

# 创建系统分类
system_groups <- list(
  "Cardiovascular" = c("Hypertension_Essential", "Hypertensive_Heart_Disease", 
                       "Angina_Pectoris", "Acute_MI", "Chronic_IHD", 
                       "Conduction_Disorders", "Atrial_Fibrillation_Flutter"),
  "Metabolic_Endocrine" = c("Type1_Diabetes", "Type2_Diabetes", "Unspecified_Diabetes",
                            "Hypothyroidism", "Hyperthyroidism", "Thyroiditis",
                            "Obesity", "Hypercholesterolemia", "Mixed_Hyperlipidemia"),
  "Gastrointestinal" = c("GERD", "Gastric_Ulcer", "Duodenal_Ulcer", "Gastritis_Duodenitis",
                         "Irritable_Bowel_Syndrome", "Constipation", 
                         "Alcoholic_Liver_Disease", "Non_alcoholic_Fatty_Liver", "Cholelithiasis"),
  "Musculoskeletal" = c("Rheumatoid_Arthritis", "Osteoarthritis", "Spondylosis",
                        "Intervertebral_Disc_Disorders", "Low_Back_Pain", "Soft_Tissue_Disorders"),
  "Mental_Behavioral" = c("Depressive_Episode", "Recurrent_Depression", 
                          "Phobic_Anxiety", "Other_Anxiety_Disorders", "Reaction_to_Severe_Stress"),
  "Respiratory" = c("Acute_URI", "Asthma", "COPD", "Pneumonia")
)

# 创建合并后的疾病分组（将低频组按系统合并）
create_consolidated_disease_groups <- function(detailed_merged_data, system_groups, threshold = 300) {
  # 识别高频疾病组（保留独立分析）
  high_freq_groups <- detailed_merged_data %>%
    filter(!Below_Threshold) %>%
    select(Disease_Group, Total_Cases)
  
  # 识别低频疾病组并按系统合并
  low_freq_by_system <- list()
  
  for (system_name in names(system_groups)) {
    system_diseases <- system_groups[[system_name]]
    
    # 找出该系统下的低频疾病
    system_low_freq <- detailed_merged_data %>%
      filter(Disease_Group %in% system_diseases & Below_Threshold) %>%
      pull(Disease_Group)
    
    if (length(system_low_freq) > 0) {
      # 获取这些低频疾病对应的ICD编码
      low_freq_icds <- c()
      for (disease in system_low_freq) {
        icds <- detailed_disease_groups[[disease]]
        low_freq_icds <- c(low_freq_icds, icds)
      }
      
      # 创建合并组
      merged_group_name <- paste0("Other_", system_name, "_Diseases")
      low_freq_by_system[[merged_group_name]] <- low_freq_icds
    }
  }
  
  # 合并高频组和新的低频合并组
  consolidated_groups <- detailed_disease_groups
  
  # 移除原来的低频组
  for (low_group in low_frequency_groups) {
    consolidated_groups[[low_group]] <- NULL
  }
  
  # 添加新的系统合并组
  for (merged_group in names(low_freq_by_system)) {
    consolidated_groups[[merged_group]] <- low_freq_by_system[[merged_group]]
  }
  
  return(consolidated_groups)
}

# 创建合并后的疾病分组
consolidated_disease_groups <- create_consolidated_disease_groups(detailed_merged_data, system_groups, 300)

# 应用合并后的分组
final_merged_data <- create_hierarchical_groups_with_threshold(icd_data, consolidated_disease_groups)

# 可视化结果
library(ggplot2)
library(scales)

# 1. 合并前后的对比
comparison_data <- data.frame(
  Stage = rep(c("Before Consolidation", "After Consolidation"), each = nrow(detailed_merged_data)),
  Disease_Group = c(detailed_merged_data$Disease_Group, final_merged_data$Disease_Group),
  Total_Cases = c(detailed_merged_data$Total_Cases, final_merged_data$Total_Cases),
  Below_Threshold = c(detailed_merged_data$Below_Threshold, final_merged_data$Below_Threshold)
)

# 只显示高频疾病组（病例数>=300）
high_freq_comparison <- comparison_data %>%
  filter(Total_Cases >= 300)

p1 <- ggplot(high_freq_comparison, aes(x = reorder(Disease_Group, Total_Cases), y = Total_Cases, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "疾病组病例数对比 (仅显示病例数≥300的组)",
    subtitle = "合并前后高频疾病组的比较",
    x = "疾病组",
    y = "病例数",
    fill = "阶段"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 2. 最终疾病组分布（按系统着色）
system_colors <- c(
  "Cardiovascular" = "#E41A1C",
  "Metabolic_Endocrine" = "#377EB8", 
  "Gastrointestinal" = "#4DAF4A",
  "Musculoskeletal" = "#984EA3",
  "Mental_Behavioral" = "#FF7F00",
  "Respiratory" = "#FFFF33"
)

# 为每个疾病组分配系统颜色
final_merged_data <- final_merged_data %>%
  mutate(
    System = case_when(
      grepl("Cardiovascular", Disease_Group) ~ "Cardiovascular",
      grepl("Metabolic|Diabetes|Thyroid|Obesity|Lipid", Disease_Group) ~ "Metabolic_Endocrine",
      grepl("Gastrointestinal|GERD|Ulcer|Bowel|Liver", Disease_Group) ~ "Gastrointestinal",
      grepl("Musculoskeletal|Arthritis|Osteoarthritis|Back", Disease_Group) ~ "Musculoskeletal",
      grepl("Mental|Depressive|Anxiety|Stress", Disease_Group) ~ "Mental_Behavioral",
      grepl("Respiratory|URI|Asthma|COPD|Pneumonia", Disease_Group) ~ "Respiratory",
      TRUE ~ "Other"
    )
  )

p2 <- ggplot(final_merged_data, aes(x = reorder(Disease_Group, Total_Cases), y = Total_Cases, fill = System)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = system_colors) +
  geom_hline(yintercept = 300, linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "最终疾病组病例数分布",
    subtitle = "红色虚线表示合并阈值 (300例)",
    x = "疾病组",
    y = "病例数",
    fill = "系统分类"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 3. 系统级别的病例数分布
system_summary <- final_merged_data %>%
  group_by(System) %>%
  summarise(
    Total_Cases = sum(Total_Cases),
    Number_of_Groups = n(),
    Mean_Cases_per_Group = round(mean(Total_Cases), 1)
  ) %>%
  arrange(desc(Total_Cases))

p3 <- ggplot(system_summary, aes(x = reorder(System, Total_Cases), y = Total_Cases, fill = System)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = system_colors) +
  geom_text(aes(label = paste0("n=", Number_of_Groups, " groups")), 
            hjust = -0.1, size = 3) +
  labs(
    title = "各系统疾病负担分布",
    x = "系统分类",
    y = "总病例数"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# 显示图表
print(p1)
print(p2)
print(p3)

# 输出统计摘要
cat("=== 合并效果统计摘要 ===\n")
cat("原始疾病组数量:", nrow(detailed_merged_data), "\n")
cat("最终疾病组数量:", nrow(final_merged_data), "\n")
cat("合并的疾病组数量:", nrow(detailed_merged_data) - nrow(final_merged_data), "\n")
cat("病例数低于300的原始组数量:", sum(detailed_merged_data$Below_Threshold), "\n")
cat("最终所有组病例数均 ≥ 300:", all(final_merged_data$Total_Cases >= 300), "\n\n")

cat("=== 最终疾病组列表 ===\n")
print(final_merged_data %>% 
        select(Disease_Group, Total_Cases, System) %>%
        arrange(desc(Total_Cases)))

# 保存最终结果
write.csv(final_merged_data, "Consolidated_Disease_Groups_Final.csv", row.names = FALSE)

# 生成用于脂肪分布分析的格式化数据
fat_analysis_ready <- final_merged_data %>%
  select(Disease_Group, Total_Cases, System) %>%
  mutate(
    Analysis_Group = ifelse(
      grepl("Other_", Disease_Group),
      paste0("Consolidated_", System),
      Disease_Group
    )
  ) %>%
  group_by(Analysis_Group, System) %>%
  summarise(
    Total_Cases = sum(Total_Cases),
    Number_of_Subgroups = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(Total_Cases))

cat("\n=== 脂肪分布分析就绪分组 ===\n")
print(fat_analysis_ready)