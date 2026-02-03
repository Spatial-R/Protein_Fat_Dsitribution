
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(ggcorrplot)
library(ggtern)
library(GGally)
library(ggVennDiagram)
library(patchwork)

source("Codes/00_Basic_Functions.R")
source("Codes_Fat/00_Fat_BasicFunctions.R")
source("Codes_Fat/00_Descritive_Functions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Fat_Dis/Pro_Impt.RData")

dat_body_fin <- filter(dat_body_fin_u,instance == 2)

tar_var <- c("Pancreas_PDFF_fat_fraction","Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration","area_of_pericardial_fat")

load("Process_Data/Fat_Dis/Model_Data.RData")
tar_col <- match(tar_var,names(dat.ana.u.1))
normal_fun <- function(data){
  (data - min(data))/(max(data) - min(data))
}
dat.ana.u.2 <- dat.ana.u.1
dat.ana.u.2[,c(tar_col)] <- apply(dat.ana.u.2[,c(tar_col)], 2, normal_fun)


main_plot <- ggtern(data = dat.ana.u.2, aes(y = Pancreas_PDFF_fat_fraction, 
                                            x = Liver_PDFF_fat_fraction, 
                                            z = Muscle_fat_infiltration)) +
  stat_density_tern(
    geom = "polygon",
    aes(fill = after_stat(level)),
    bins = 25,color = NA,
    alpha = 0.9,bdl = 0.05,linetype = "dashed",
    bdl.val = 0.005    # 将低于检测限的值设置为0.005（在密度计算中处理）
  ) +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  #annotate("text", x = 0.5, y = -0.4, z = 0.1, label = "Custom X Label", color = "blue") +
  #annotate("text", x = 0.1, y = 0.5, z = 0.4, label = "Custom Y Label", color = "darkgreen") +
  #annotate("text", x = 0.3, y = 0.1, z = 0.6, label = "Custom Z Label", color = "purple")+
  geom_point(aes(color = area_of_pericardial_fat),alpha = 0.8,size = 0.5) +
  geom_Tline(Tintercept = 5, 
             color = "red", alpha = 1, size = 0.5) +
  geom_Lline(Lintercept = 8.3, 
             color = "blue", alpha = 1, size = 0.5) +
  geom_Rline(Rintercept = 9, 
             color = "green", alpha = 1, size = 0.5) +
  scale_color_viridis_c(
    option = "viridis",
    name = "Area of pericardial fat (cm²)",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(3, "cm"),
      barheight = unit(0.3, "cm"),
      order = 1 
    )
  ) +
  labs(
   # title = "Multivariate Analysis of Ectopic Fat Deposition",
    #subtitle = "Relationship between Pancreatic, Hepatic, and Muscular Fat Content",
    x = "Liver \nPDFF(%)", 
    y = "Pancreas \nPDFF(%)",
    z = "Muscle \nFI(%)") +
  
  theme_rgbg(base_size = 10, base_family = "serif") +
  theme(
    legend.position = "bottom",
    plot.margin = unit(rep(2,4),"cm"),
    legend.box = "vertical",  # 垂直排列图例
    legend.spacing = unit(0.2, "cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    tern.panel.background = element_rect(fill = "white", color = "black"), 
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    plot.caption = element_text(size = 9, color = "gray50"),
    axis.title = element_text(face = "bold"))

ggsave(filename = "Figures/Fat/Basic_Vis/Basic_ggtern.tiff", main_plot,
       width = 20,height = 18,units = "cm",dpi = 300)

marginal_plots <- create_marginal_plots(dat.ana.u.1,size_u = 12)

marginal_fig <- cowplot::plot_grid(plotlist = marginal_plots,labels = c(""),ncol = 4)
# 
# marginal_composite <- patchwork::wrap_plots(marginal_plots, ncol = 1) +
#   plot_annotation(title = "Marginal Distributions of Fat Deposition Variables",
#                   theme = theme(plot.title = element_text(hjust = 0.5, size = 12)))


##########################  Sex and age distribution ###########################

fat_demo_dat <- dat.ana.u.1[, match(c(tar_var,"age","sex"),names(dat.ana.u.1))]
names(fat_demo_dat)[1:4] <- fat_abbr

fat_demo_dat <- mutate(fat_demo_dat,
                      sex = ifelse(sex == 1,"Male","Female"),
                      Age_Group = ifelse(age >= 60, "≥60 years", "<60 years"),
                      Age_Group = factor(Age_Group, levels = c("<60 years", "≥60 years"))) # Set factor order)

fat_demo_dat_u <- fat_demo_dat[,-c(5)] %>%
  pivot_longer(cols = c(fat_abbr),names_to = "Organ",values_to = "Fat_Percentage")

fat_demo_dat_long <- fat_demo_dat[,-c(5)] %>%
  pivot_longer(
    cols = all_of(fat_abbr),
    names_to = "Organ",
    values_to = "Fat_Percentage"
  ) %>%
  mutate(
    Organ = factor(Organ, levels = fat_abbr)
  )

# 4. Create panels for each organ
organ_panels <- list()

for (organ in fat_abbr) {
  # Create left panel (age comparison)
  age_panel <- create_age_panel(organ, fat_demo_dat_long)
  
  # Create right panel (sex comparison)
  sex_panel <- create_sex_panel(organ, fat_demo_dat_long)
  
  # Combine into a list
  organ_panels[[paste0(organ, "_age")]] <- age_panel
  organ_panels[[paste0(organ, "_sex")]] <- sex_panel
}

# First, create a list of panel pairs (each organ: age + sex)
panel_rows <- list()

for (organ in fat_abbr) {
  # Get the two panels for this organ
  age_panel <- organ_panels[[paste0(organ, "_age")]]
  sex_panel <- organ_panels[[paste0(organ, "_sex")]]
  
  # Combine them into one row
  combined_row <- age_panel + sex_panel
  
  # Store the combined row
  panel_rows[[organ]] <- combined_row
}

# 6. Create combined plot with compressed spacing
final_plot <- wrap_plots(panel_rows, ncol = 1, heights = rep(1, length(fat_abbr))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")

fig_a <- add_column_labels(final_plot)


age_comparison_stats <- data.frame()
sex_comparison_stats <- data.frame()

for (organ in fat_abbr) {
  # Filter data for current organ
  organ_data <- fat_demo_dat_long %>% filter(Organ == organ)
  
  # Age comparison
  age_test <- t.test(Fat_Percentage ~ Age_Group, data = organ_data)
  age_stats <- organ_data %>%
    group_by(Age_Group) %>%
    summarise(
      N = n(),
      Mean = mean(Fat_Percentage, na.rm = TRUE),
      SD = sd(Fat_Percentage, na.rm = TRUE),
      .groups = "drop"
    )
  
  age_comparison_stats <- rbind(
    age_comparison_stats,
    data.frame(
      Organ = organ,
      Group1 = "<60 years",
      Group2 = "≥60 years",
      N1 = age_stats$N[1],
      N2 = age_stats$N[2],
      Mean1 = round(age_stats$Mean[1], 2),
      Mean2 = round(age_stats$Mean[2], 2),
      SD1 = round(age_stats$SD[1], 2),
      SD2 = round(age_stats$SD[2], 2),
      p_value = round(age_test$p.value, 4)
    )
  )
  
  # Sex comparison
  sex_test <- t.test(Fat_Percentage ~ sex, data = organ_data)
  sex_stats <- organ_data %>%
    group_by(sex) %>%
    summarise(
      N = n(),
      Mean = mean(Fat_Percentage, na.rm = TRUE),
      SD = sd(Fat_Percentage, na.rm = TRUE),
      .groups = "drop"
    )
  
  sex_comparison_stats <- rbind(
    sex_comparison_stats,
    data.frame(
      Organ = organ,
      Group1 = "Male",
      Group2 = "Female",
      N1 = sex_stats$N[sex_stats$sex == "Male"],
      N2 = sex_stats$N[sex_stats$sex == "Female"],
      Mean1 = round(sex_stats$Mean[sex_stats$sex == "Male"], 2),
      Mean2 = round(sex_stats$Mean[sex_stats$sex == "Female"], 2),
      SD1 = round(sex_stats$SD[sex_stats$sex == "Male"], 2),
      SD2 = round(sex_stats$SD[sex_stats$sex == "Female"], 2),
      p_value = round(sex_test$p.value, 4)
    )
  )
}

cat("\n=== AGE COMPARISON STATISTICS ===\n")
print(age_comparison_stats, row.names = FALSE)

cat("\n=== SEX COMPARISON STATISTICS ===\n")
print(sex_comparison_stats, row.names = FALSE)

write.csv(age_comparison_stats, "Figures/Fat/Basic_Vis/Age_Comparison_Stats.csv", row.names = FALSE)
write.csv(sex_comparison_stats, "Figures/Fat/Basic_Vis/Sex_Comparison_Stats.csv", row.names = FALSE)


# load(file = "Figures/Fat/Clust.RData")
# names(dat_mclust_tem) <- gsub("_residual","",names(dat_mclust_tem))
# names(dat_mclust_tem)[1:4] <-  fat_abbr
# dat_mclust_tem <- mutate(dat_mclust_tem,
#                          cluster = factor(cluster,levels = c(1:4),
#                                           labels = paste0("Cluster ",1:4)))

matrix_dat <- dat.ana.u.1[, match(c(tar_var,"sex"),names(dat.ana.u.1))]
matrix_dat <- mutate(matrix_dat,sex = ifelse(sex == 1,"Male","Female"))
names(matrix_dat)[1:4] <- fat_abbr
fig_b <- ggpairs(matrix_dat,
                 columns = c(fat_abbr),
                 mapping = ggplot2::aes(color = sex, alpha = 0.5),
                 upper = list(continuous = wrap("cor", size = 4)),
                 lower = list(continuous = wrap("points", size = 1.5)),
                 diag = list(continuous = wrap("densityDiag", alpha = 0.7))) +
  theme_bw(base_size = 12,base_family = "serif")

ggsave(filename = "Figures/Fat/Basic_Vis/Basic_ggparis.tiff",fig_b,width = 20,
       height = 16,units = "cm",dpi = 300)


g1 <- GGally::ggmatrix_gtable(fig_b)
fig_merge_1 <- plot_grid(fig_a,g1,ncol = 2,labels = c("A","B"))


library(psych)

cor_dat <- fat_demo_dat[,1:4]
names(cor_dat) <- c("Pancrea", "Liver", "Muscle", "Card")
cor_results <- corr.test(cor_dat, use = "pairwise", method = "pearson", adjust = "none")
cor_results$ci

correlation_plot <- create_correlation_heatmap(dat.ana.u.1)


cardiac_thresholds <- dat.ana.u.1 %>%
  group_by(sex) %>%
  summarise(
    thr_car = quantile(area_of_pericardial_fat, 0.75, na.rm = TRUE),
    muc_car = quantile(Muscle_fat_infiltration, 0.75, na.rm = TRUE), .groups = 'drop')

dat.ana.u.2 <- dat.ana.u.1 %>%
  left_join(cardiac_thresholds, by = "sex") %>%
  mutate(cardiac_fat_abnormal = area_of_pericardial_fat > thr_car,
         pancreatic_fat = Pancreas_PDFF_fat_fraction >= log(10.4),
         hepatic_steatosis = Liver_PDFF_fat_fraction > log(5.0),
         muscular_fat_infiltration = Muscle_fat_infiltration > muc_car)

sets <- list(
  Pericardial = dat.ana.u.2$eid[as.integer(dat.ana.u.2$cardiac_fat_abnormal) == 1],
  Liver = dat.ana.u.2$eid[as.integer(dat.ana.u.2$hepatic_steatosis) == 1],
  Pancreas = dat.ana.u.2$eid[as.integer(dat.ana.u.2$pancreatic_fat) == 1],
  Muscle = dat.ana.u.2$eid[as.integer(dat.ana.u.2$muscular_fat_infiltration) == 1])

fig_3 <- ggVennDiagram(sets, 
                       category.names = c("PFA","Liver PDFF", "Pancreas PDFF", 
                                          "Muscle FI"),
                       label = "both",  # Display both count and percentage
                       label_size = 3.5,label_alpha = 0,set_size = 4,
                       edge_size = 1.2) +
  # Option 1: Use viridis with lighter range
  scale_fill_viridis_c(option = "plasma", 
                       direction = -1,
                       begin = 0.3,    # Start with lighter color
                       end = 0.8,      # End with medium depth
                       alpha = 0.7,    # Increase transparency
                       name = "Number of Overlapping Samples",
                       guide = guide_colorbar(
                         barwidth = 12,
                         barheight = 0.8,
                         title.position = "top",
                         title.hjust = 0.5)) +
  
  # Option 2: Use softer color scheme
  # scale_fill_gradient(low = "#E6F2FF", high = "#1E88E5",  # Light blue to medium blue
  #                     name = "Number of Overlapping Samples") +
  
  # Option 3: Use brewer palette
  # scale_fill_distiller(palette = "Blues", direction = 1,
  #                      name = "Number of Overlapping Samples") +
  
  # Use lighter border colors
  scale_color_manual(values = c("#B0BEC5", "#B0BEC5", "#B0BEC5", "#B0BEC5")) +  # Light gray
  
  # Professional theme
  theme_void(base_family = "serif") +
  theme(
    plot.title = element_text(
      hjust = 0.5,size = 16, face = "bold",color = "#4A6572"),  # Lighter title color
    plot.subtitle = element_text(
      hjust = 0.5,size = 10,color = "#78909C"),  # Lighter subtitle color
    plot.caption = element_text(
      hjust = 0.5,size = 10,color = "#90A4AE"),  # Lighter caption color
    legend.position = "bottom",
    legend.title = element_text(
      size = 10,face = "bold",color = "#546E7A"),  # Lighter legend title color
    legend.text = element_text(
      size = 8,color = "#607D8B"),  # Lighter legend text color
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Optional: Add light background to plot area
    # panel.background = element_rect(fill = "#F5F7FA", color = NA),
    # Optional: Adjust legend background
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA))

fig_r <- plot_grid(correlation_plot,
                   fig_3,ncol = 1,labels = c("B","C","D"))

fig_merge <- plot_grid(fig_a,fig_r,ncol = 2,labels = c("A",""))

ggsave(filename = "Figures/Fat/Basic_Vis/Basic_Vis.tiff",fig_merge,
       width = 28,height = 20,units = "cm",dpi = 300)





total_participants <- nrow(dat.ana.u.2)

pericardial_only <- sets$Pericardial[!sets$Pericardial %in% 
                                       c(sets$Liver, sets$Pancreas, sets$Muscle)]
liver_only <- sets$Liver[!sets$Liver %in% 
                           c(sets$Pericardial, sets$Pancreas, sets$Muscle)]
pancreas_only <- sets$Pancreas[!sets$Pancreas %in% 
                                 c(sets$Pericardial, sets$Liver, sets$Muscle)]
muscle_only <- sets$Muscle[!sets$Muscle %in% 
                             c(sets$Pericardial, sets$Liver, sets$Pancreas)]

liver_pancreas <- intersect(sets$Liver, sets$Pancreas)
liver_muscle <- intersect(sets$Liver, sets$Muscle)
liver_pericardial <- intersect(sets$Liver, sets$Pericardial)
pancreas_muscle <- intersect(sets$Pancreas, sets$Muscle)
pancreas_pericardial <- intersect(sets$Pancreas, sets$Pericardial)
muscle_pericardial <- intersect(sets$Muscle, sets$Pericardial)

liver_pancreas_only <- liver_pancreas[!liver_pancreas %in% 
                                        c(sets$Muscle, sets$Pericardial)]
liver_muscle_only <- liver_muscle[!liver_muscle %in% 
                                    c(sets$Pancreas, sets$Pericardial)]
liver_pericardial_only <- liver_pericardial[!liver_pericardial %in% 
                                              c(sets$Pancreas, sets$Muscle)]
pancreas_muscle_only <- pancreas_muscle[!pancreas_muscle %in% 
                                          c(sets$Liver, sets$Pericardial)]
pancreas_pericardial_only <- pancreas_pericardial[!pancreas_pericardial %in% 
                                                    c(sets$Liver, sets$Muscle)]
muscle_pericardial_only <- muscle_pericardial[!muscle_pericardial %in% 
                                                c(sets$Liver, sets$Pancreas)]

liver_pancreas_muscle <- Reduce(intersect, list(sets$Liver, sets$Pancreas, sets$Muscle))
liver_pancreas_pericardial <- Reduce(intersect, list(sets$Liver, sets$Pancreas, sets$Pericardial))
liver_muscle_pericardial <- Reduce(intersect, list(sets$Liver, sets$Muscle, sets$Pericardial))
pancreas_muscle_pericardial <- Reduce(intersect, list(sets$Pancreas, sets$Muscle, sets$Pericardial))

liver_pancreas_muscle_only <- liver_pancreas_muscle[!liver_pancreas_muscle %in% 
                                                      sets$Pericardial]
liver_pancreas_pericardial_only <- liver_pancreas_pericardial[!liver_pancreas_pericardial %in% 
                                                                sets$Muscle]
liver_muscle_pericardial_only <- liver_muscle_pericardial[!liver_muscle_pericardial %in% 
                                                            sets$Pancreas]
pancreas_muscle_pericardial_only <- pancreas_muscle_pericardial[!pancreas_muscle_pericardial %in% 
                                                                  sets$Liver]

all_four <- Reduce(intersect, list(sets$Pericardial, sets$Liver, 
                                   sets$Pancreas, sets$Muscle))

any_abnormality <- unique(c(sets$Pericardial, sets$Liver, sets$Pancreas, sets$Muscle))
any_abnormality_count <- length(any_abnormality)

count_compartments <- function(eid) {
  sum(eid %in% sets$Pericardial,
      eid %in% sets$Liver,
      eid %in% sets$Pancreas,
      eid %in% sets$Muscle)
}

compartment_counts <- sapply(any_abnormality, count_compartments)
multiple_compartments_count <- sum(compartment_counts >= 2)

venn_combination_table <- data.frame(
  Combination = c(
    "Pericardial only",
    "Liver only",
    "Pancreas only",
    "Muscle only",
    
    "Liver & Pancreas only",
    "Liver & Muscle only",
    "Liver & Pericardial only",
    "Pancreas & Muscle only",
    "Pancreas & Pericardial only",
    "Muscle & Pericardial only",
    
    "Liver, Pancreas & Muscle only",
    "Liver, Pancreas & Pericardial only",
    "Liver, Muscle & Pericardial only",
    "Pancreas, Muscle & Pericardial only",
    
    "All four compartments",
    
    "Any abnormality (total)",
    "Multiple compartments (≥2)",
    "Total participants"),
  Count = c(
    length(pericardial_only),
    length(liver_only),
    length(pancreas_only),
    length(muscle_only),
    
    length(liver_pancreas_only),
    length(liver_muscle_only),
    length(liver_pericardial_only),
    length(pancreas_muscle_only),
    length(pancreas_pericardial_only),
    length(muscle_pericardial_only),
    
    length(liver_pancreas_muscle_only),
    length(liver_pancreas_pericardial_only),
    length(liver_muscle_pericardial_only),
    length(pancreas_muscle_pericardial_only),
    length(all_four),
    any_abnormality_count,
    multiple_compartments_count,
    total_participants),
  
  Percentage = c(
    round(length(pericardial_only) / total_participants * 100, 2),
    round(length(liver_only) / total_participants * 100, 2),
    round(length(pancreas_only) / total_participants * 100, 2),
    round(length(muscle_only) / total_participants * 100, 2),
    
    round(length(liver_pancreas_only) / total_participants * 100, 2),
    round(length(liver_muscle_only) / total_participants * 100, 2),
    round(length(liver_pericardial_only) / total_participants * 100, 2),
    round(length(pancreas_muscle_only) / total_participants * 100, 2),
    round(length(pancreas_pericardial_only) / total_participants * 100, 2),
    round(length(muscle_pericardial_only) / total_participants * 100, 2),
    
    round(length(liver_pancreas_muscle_only) / total_participants * 100, 2),
    round(length(liver_pancreas_pericardial_only) / total_participants * 100, 2),
    round(length(liver_muscle_pericardial_only) / total_participants * 100, 2),
    round(length(pancreas_muscle_pericardial_only) / total_participants * 100, 2),
    
    round(length(all_four) / total_participants * 100, 2),
    
    round(any_abnormality_count / total_participants * 100, 2),
    round(multiple_compartments_count / total_participants * 100, 2),
    100.00
  ),
  stringsAsFactors = FALSE
)

venn_combination_table$Category <- c(
  rep("Single compartment", 4),
  rep("Two compartments", 6),
  rep("Three compartments", 4),
  "Four compartments",
  rep("Summary statistics", 3)
)

venn_combination_table <- venn_combination_table %>%
  select(Category, Combination, Count, Percentage)

cat("\n=== COMPLETE VENN COMBINATION TABLE ===\n")
print(venn_combination_table, row.names = FALSE)

write.csv(venn_combination_table, 
          "Figures/Fat/Basic_Vis/Venn_Combination_Table.csv", 
          row.names = FALSE)

non_zero_table <- venn_combination_table %>%
  filter(Count > 0)

cat("\n=== NON-ZERO COMBINATIONS ONLY ===\n")
print(non_zero_table, row.names = FALSE)

category_summary <- venn_combination_table %>%
  filter(Category != "Summary statistics") %>%
  group_by(Category) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Percentage = sum(Percentage),
    Number_of_Combinations = n()
  )

cat("\n=== CATEGORY SUMMARY ===\n")
print(category_summary, row.names = FALSE)

sum_of_combinations <- sum(venn_combination_table$Count[venn_combination_table$Category != "Summary statistics"])
sum_of_summary <- venn_combination_table$Count[venn_combination_table$Combination == "Any abnormality (total)"]

cat("\n=== DATA VALIDATION ===\n")
cat(sprintf("Sum of all combinations: %d\n", sum_of_combinations))
cat(sprintf("Total with any abnormality: %d\n", sum_of_summary))
cat(sprintf("Difference: %d\n", sum_of_combinations - sum_of_summary))

if (abs(sum_of_combinations - sum_of_summary) < 1) {
  cat("\nData validation PASSED: Sum of all combinations equals total with any abnormality.\n")
  
  cat("\n=== SAMPLE RESULT DESCRIPTION ===\n")
  cat(sprintf("The Venn analysis of four body compartments revealed that %d participants (%.1f%%) had at least one abnormal fat deposition pattern. ", 
              any_abnormality_count, 
              round(any_abnormality_count/total_participants*100, 1)))
  cat(sprintf("Among these, %d participants (%.1f%%) exhibited abnormalities in multiple compartments. ", 
              multiple_compartments_count, 
              round(multiple_compartments_count/total_participants*100, 1)))
  
  top_combinations <- non_zero_table %>%
    filter(!Combination %in% c("Any abnormality (total)", "Multiple compartments (≥2)", "Total participants")) %>%
    arrange(desc(Count))
  
  if (nrow(top_combinations) > 0) {
    cat(sprintf("The most frequent pattern was '%s' (n=%d, %.1f%%), ", 
                top_combinations$Combination[1], 
                top_combinations$Count[1],
                top_combinations$Percentage[1]))
    
    if (nrow(top_combinations) > 1) {
      cat(sprintf("followed by '%s' (n=%d, %.1f%%). ", 
                  top_combinations$Combination[2], 
                  top_combinations$Count[2],
                  top_combinations$Percentage[2]))
    } else {
      cat(". ")
    }
  }
  cat("These patterns suggest a systemic nature of abnormal fat distribution in this population.\n")
} else {
  cat("\nData validation FAILED: Please check the calculations.\n")
}


plot_data <- category_summary %>%
  filter(Category != "Summary statistics")

plot_data <- mutate(plot_data,Category = factor(Category,
                                                levels = c("Single compartment",
                                                           "Two compartments",
                                                           "Three compartments",
                                                           "Four compartments")))

fig_summ <- ggplot(plot_data, aes(x = Category, y = Total_Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", Total_Percentage)), 
            vjust = -0.5, size = 4) +
  labs(title = "Distribution of Abnormal Fat Deposition Patterns",
       x = "Number of Affected Compartments",
       y = "Percentage of Participants (%)") +
  theme_classic(base_size = 12,base_family = "serif") +
  theme(legend.position = "none")

ggsave("Figures/Fat/Basic_Vis/Venn_Category_Summary.png", 
       fig_summ, width = 8, height = 6, dpi = 300)
