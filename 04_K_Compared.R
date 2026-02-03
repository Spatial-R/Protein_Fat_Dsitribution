
rm(list = ls())

library(dplyr)
library(ggVennDiagram)
library(ggplot2)
library(gridExtra)
library(grid)

protein_fin <- data.frame(); fig_list <- list()

for (tar_id in (1:4)){

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")
dir_path <- "Process_Data/Fat_Dis/Fat_Causes_Res_NM/"

load(paste0(dir_path,tar_var[tar_id],"_7_100_kocmi_nm.RData"))
tar_res_cause_7 <- filter(res_kocmi,target %in% "target" & p_adj < 0.05)
tar_res_out_7 <- filter(res_kocmi,regulator %in% "target" & p_adj < 0.05)
print(tar_res_cause_7);print(tar_res_out_7);
tar_res_cause_7 <- mutate(tar_res_cause_7,id = tar_var[tar_id])
protein_fin <- rbind(protein_fin,tar_res_cause_7)

load(paste0(dir_path,tar_var[tar_id],"_3_100_kocmi_nm.RData"))
tar_res_cause_3 <- filter(res_kocmi,target %in% "target" & p_adj < 0.05)
tar_res_out_3 <- filter(res_kocmi,regulator %in% "target" & p_adj < 0.05)
print(tar_res_cause_3);print(tar_res_out_3);

load(paste0(dir_path,tar_var[tar_id],"_5_100_kocmi_nm.RData"))
tar_res_cause_5 <- filter(res_kocmi,target %in% "target" & p_adj < 0.05)
tar_res_out_5 <- filter(res_kocmi,regulator %in% "target" & p_adj < 0.05)
print(tar_res_cause_5);print(tar_res_out_5);


# Create data
set1 <- toupper(tar_res_cause_3$regulator)
set2 <- toupper(tar_res_cause_5$regulator)
set3 <- toupper(tar_res_cause_7$regulator)

my_sets <- list(Dataset_A = set1, Dataset_B = set2, Dataset_C = set3)

# Function to calculate Venn diagram elements
calculate_venn_elements <- function(A, B, C) {
  list(
    only_A = setdiff(A, union(B, C)),
    only_B = setdiff(B, union(A, C)),
    only_C = setdiff(C, union(A, B)),
    AB_only = setdiff(intersect(A, B), C),
    AC_only = setdiff(intersect(A, C), B),
    BC_only = setdiff(intersect(B, C), A),
    ABC = Reduce(intersect, list(A, B, C))
  )
}

# Calculate elements for each region
venn_elements <- calculate_venn_elements(set1, set2, set3)

# Create Venn diagram
venn_plot <- ggVennDiagram(my_sets, label = "count",
                           category.names = c("k = 3", "k = 5", "k = 7")) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  theme(legend.position = "none") #+
  #labs(title = "Protein Distribution Under Three k-Value Settings (3, 5, 7)")

# Create element table with line breaks every 6 proteins
element_table <- data.frame(
  Region = c("K3-specific", "K5-specific", "K7-specific",
             "K3-K5 Shared (K7-exclusive)", "K3-K7 Shared (K5-exclusive)", 
             "K5-K7 Shared (K3-exclusive)", "Core Intersection (All Conditions)"),
  Proteins = sapply(venn_elements, function(x) {
    # Add line break every 6 proteins
    if(length(x) > 0) {
      # Split proteins into groups of 6
      groups <- split(x, ceiling(seq_along(x)/6))
      # Combine each group with line breaks
      paste(sapply(groups, function(g) paste(g, collapse = ", ")), collapse = ",\n")
    } else {
      ""
    }
  }),
  Count = sapply(venn_elements, length)
)
row.names(element_table) <- NULL

# Print table
print(element_table)

# Create combined plot with adjusted table
combined_plot <- grid.arrange(
  venn_plot,
  tableGrob(element_table, rows = NULL, 
            theme = ttheme_default(base_size = 8)),
  nrow = 2,
  heights = c(1, 1)
)
fig_list[[tar_id]] <- combined_plot
ggsave(paste0("Figures/Fat/K_Sensitive/",tar_var[tar_id],"_venn.png"), 
       plot = combined_plot,
       width = 14, height = 10,dpi = 300)
}

library(cowplot)

# 将grob对象转换为ggplot对象
fig_list_ggplot <- lapply(fig_list, function(x) as_ggplot(x))

# 使用plot_grid合并，并添加标签
fig_merge <- plot_grid(plotlist = fig_list_ggplot, 
                       nrow = 2, ncol = 2, 
                       labels = c("A", "B", "C", "D"), 
                       label_size = 12)
ggsave("Figures/Fat/K_Sensitive/merged_venn.png", 
       plot = fig_merge,
       width = 14, height = 12, dpi = 300)
