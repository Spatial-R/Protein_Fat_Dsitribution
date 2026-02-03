
create_marginal_plots <- function(data,size_u = 12) {
  
  p1 <- ggplot(data, aes(x = Liver_PDFF_fat_fraction)) +
    geom_density(fill = "#08306B", alpha = 0.7) +
    labs(x = "Liver PDFF (%)", y = "Density") +
    theme_minimal(base_size = size_u) +
    theme(axis.title = element_text(size = size_u))
  
  p2 <- ggplot(data, aes(x = Pancreas_PDFF_fat_fraction)) +
    geom_density(fill = "#2171B5", alpha = 0.7) +
    labs(x = "Pancreas PDFF (%)", y = "Density") +
    theme_minimal(base_size = size_u) +
    theme(axis.title = element_text(size = size_u))
  
  p3 <- ggplot(data, aes(x = Muscle_fat_infiltration)) +
    geom_density(fill = "#4292C6", alpha = 0.7) +
    labs(x = "Muscle FI (%)", y = "Density") +
    theme_minimal(base_size = size_u) +
    theme(axis.title = element_text(size = size_u))
  
  p4 <- ggplot(data, aes(x = area_of_pericardial_fat)) +
    geom_density(fill = "#6BAED6", alpha = 0.7) +
    labs(x = "Area of pericardial fat (cm2)", y = "Density") +
    theme_minimal(base_size = size_u) +
    theme(axis.title = element_text(size = size_u))
  
  return(list(p1, p2, p3, p4))
}


create_age_panel <- function(organ_name, data) {
  # Filter data for specific organ
  organ_data <- data %>% filter(Organ == organ_name)

  age_test <- t.test(Fat_Percentage ~ Age_Group, data = organ_data)
  
  # Format p-value
  format_p_value <- function(p) {
    if (p < 0.001) return("<0.001")
    if (p < 0.01) return(sprintf("%.3f", p))
    if (p < 0.05) return(sprintf("%.3f", p))
    return(sprintf("%.3f", p))
  }
  
  age_p_label <- format_p_value(age_test$p.value)
  
  # Create age comparison plot
  p <- ggplot(organ_data, aes(x = Age_Group, y = Fat_Percentage, fill = Age_Group)) +
    # Boxplot
    geom_boxplot(
      alpha = 0.8,
      outlier.shape = 16,
      outlier.size = 1,
      outlier.alpha = 0.3,
      width = 0.6
    ) +
    # Points for distribution
    geom_jitter(
      position = position_jitter(width = 0.2),
      alpha = 0.2,
      size = 0.5,
      color = "gray40"
    ) +
    # Color settings
    scale_fill_manual(
      values = c("<60 years" = "#4E79A7", "≥60 years" = "#E15759"),
      name = "Age Group"
    ) +
    # Add significance annotation
    annotate(
      "text",
      x = 1.5,
      y = max(organ_data$Fat_Percentage, na.rm = TRUE) * 0.85,
      label = paste("p: ",age_p_label),
      size = 4,
      fontface = "bold",
      color = "darkblue"
    ) +
    labs(
      x = NULL,  
      y = organ_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5)#,
      #plot.margin = margin(c(2, 5, 2, 5), unit = "pt")  
    )
  
  return(p)
}

# 3. Create function to generate right panel (sex comparison)
create_sex_panel <- function(organ_name, data) {
  # Filter data for specific organ
  organ_data <- data %>% filter(Organ == organ_name)
  
  # Sex comparison (ignoring age)
  sex_test <- t.test(Fat_Percentage ~ sex, data = organ_data)
  
  # Format p-value
  format_p_value <- function(p) {
    if (p < 0.001) return("<0.001")
    if (p < 0.01) return(sprintf("%.3f", p))
    if (p < 0.05) return(sprintf("%.3f", p))
    return(sprintf("%.3f", p))
  }
  
  sex_p_label <- format_p_value(sex_test$p.value)
  
  # Create sex comparison plot
  p <- ggplot(organ_data, aes(x = sex, y = Fat_Percentage, fill = sex)) +
    # Boxplot
    geom_boxplot(
      alpha = 0.8,
      outlier.shape = 16,
      outlier.size = 1,
      outlier.alpha = 0.3,
      width = 0.6
    ) +
    # Points for distribution
    geom_jitter(
      position = position_jitter(width = 0.2),
      alpha = 0.2,
      size = 0.5,
      color = "gray40"
    ) +
    # Color settings
    scale_fill_manual(
      values = c("Male" = "#4E79A7", "Female" = "#E15759"),
      name = "Sex"
    ) +
    # Add significance annotation
    annotate(
      "text",
      x = 1.5,
      y = max(organ_data$Fat_Percentage, na.rm = TRUE) * 0.85,
      label = paste("p: ",sex_p_label),
      size = 4,
      fontface = "bold",
      color = "darkred"
    ) +
   
    labs(
      x = NULL,  
      y = NULL   
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5)#,
      #plot.margin = margin(c(2, 5, 2, 5), unit = "pt")  
    )
  
  return(p)
}


add_column_labels <- function(plot_obj) {

  age_label_plot <- ggplot() + 
    theme_void() +
    labs(x = "Age Group") +
    theme(
      axis.title.x = element_text(size = 12, face = "bold", vjust = 1, hjust = 0.5)
    )
  
  sex_label_plot <- ggplot() + 
    theme_void() +
    labs(x = "Sex") +
    theme(
      axis.title.x = element_text(size = 12, face = "bold", vjust = 1, hjust = 0.5)
    )
  
  label_row <- age_label_plot + sex_label_plot

  combined_plot <- plot_grid(plot_obj,label_row,ncol = 1,rel_heights = c(0.95,0.05))
  return(combined_plot)
}


create_correlation_heatmap <- function(data) {
  
  cor_vars <- data[, match(tar_var,names(data))]
  colnames(cor_vars) <- fat_abbr
  
  cor_matrix <- cor(cor_vars, use = "complete.obs")
  melted_cor <- reshape2::melt(cor_matrix)
  
  p <- ggplot(melted_cor, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), 
              color = "black", size = 3.5, fontface = "bold") +
    scale_fill_gradient2(
      high = "#D73027", low = "#4575B4", mid = "white",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = "Correlation"
    ) +
    labs(x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-10, 0, 0, 0),  # 减少图例上方边距
      plot.margin = margin(10, 10, 5, 10)  # 减少底部边距
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(4, "cm"),
        barheight = unit(0.3, "cm"),
        direction = "horizontal"
      )
    )
  return(p)
}
