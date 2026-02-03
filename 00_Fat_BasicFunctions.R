parse_icd_range <- function(s) {
  s <- stri_trim_both(s)
  if (stri_detect_fixed(s, "-")) {
    parts <- stri_split_fixed(s, "-")[[1]]
    start <- parts[1]
    end <- parts[2]
    
    letter_start <- stri_extract_first_regex(start, "^[A-Za-z]+")
    letter_end <- stri_extract_first_regex(end, "^[A-Za-z]+")
    
    if (letter_start != letter_end) {
      warning("Different letter parts in range: ", s)
    }
    
    num_str_start <- stri_extract_first_regex(start, "\\d+")
    num_str_end <- stri_extract_first_regex(end, "\\d+")
    
    n_start <- nchar(num_str_start)
    n_end <- nchar(num_str_end)
    
    num_start <- as.integer(num_str_start)
    num_end <- as.integer(num_str_end)
    
    nums <- num_start:num_end
    
    if (n_start == n_end) {
      formatted_nums <- stri_pad_left(nums, width = n_start, pad = "0")
      return(paste0(letter_start, formatted_nums))
    } else {
      warning("Different digit lengths in range: ", s, ". Using integer representation without padding.")
      return(paste0(letter_start, nums))
    }
  } else {
    return(s)
  }
}


# Function to create binary outcomes based on gender-specific thresholds
create_binary_outcomes <- function(data, thresholds) {
  
  # Ensure sex is coded properly (1 = Male, 0 = Female)
  if(!all(unique(data$sex) %in% c("0", "1"))) {
    stop("Sex variable should be coded as 1 for Male and 0 for Female")
  }
  
  # Create binary outcomes for each fat measurement
  data <- data %>%
    mutate(
      # Pancreas PDFF
      Pancreas_PDFF_binary = as.factor(
        ifelse(sex == 1,  # Male
               ifelse(Pancreas_PDFF_fat_fraction > clinical_thresholds$Male$Pancreas_PDFF, 1, 0),
               ifelse(Pancreas_PDFF_fat_fraction > clinical_thresholds$Female$Pancreas_PDFF, 1, 0)
        )
      ),
      
      # Liver PDFF
      Liver_PDFF_binary = as.factor(
        ifelse(sex == 1,  # Male
               ifelse(Liver_PDFF_fat_fraction > clinical_thresholds$Male$Liver_PDFF, 1, 0),
               ifelse(Liver_PDFF_fat_fraction > clinical_thresholds$Female$Liver_PDFF, 1, 0)
        )
      ),
      
      # Muscle fat infiltration
      Muscle_fat_binary = as.factor(
        ifelse(sex == 1,  # Male
               ifelse(Muscle_fat_infiltration > clinical_thresholds$Male$Muscle, 1, 0),
               ifelse(Muscle_fat_infiltration > clinical_thresholds$Female$Muscle, 1, 0)
        )
      ),
      
      # Pericardial fat area
      Pericardial_fat_binary = as.factor(
        ifelse(sex == 1,  # Male
               ifelse(area_of_pericardial_fat > clinical_thresholds$Male$area_of_peri, 1, 0),
               ifelse(area_of_pericardial_fat > clinical_thresholds$Female$area_of_peri, 1, 0)
        )
      )
    )
  return(data)
}


fast_impute <- function(data, method = "median", k = 5) {
  require(dplyr)
  if (method == "median") {
    data <- data %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
  } else if (method == "mean") {
    data <- data %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  } else if (method == "knn") {
    require(VIM)
    data <- kNN(data, k = k)[, 1:ncol(data)]
  } else if (method == "rf") {
    require(missRanger)
    data <- missRanger(data, num.trees = 50, maxiter = 3)
  }
  return(data)
}


fat_abbr <- c("Pancreas PDFF","Liver PDFF","Muscle FI", "PFA")


batch_remove_effects <- function(data, outcome_vars, 
                                 covariate_vars = c("age", "sex"),
                                 nonliner_vars = c("age"),
                                 model_type = "lm") {

  residuals_list <- lapply(outcome_vars, function(var) {
    data$outcome <- data[[var]]
    if(model_type == "lm"){
      formula_str <- paste("outcome ~", paste(covariate_vars, collapse = " + "))
      model <- lm(as.formula(formula_str), data = data) 
    } else {
      non_form <- unlist(lapply(nonliner_vars,function(idd){
        paste0("s(",idd,")")
      }))
      vars_tem <- covariate_vars[-c(match(nonliner_vars,covariate_vars))]
      formula_str <- paste("outcome ~", paste(c(vars_tem,non_form), collapse = " + "))
      model <- mgcv::gam(as.formula(formula_str), data = data)  
    }
    residuals(model)
  })
  residuals_df <- as.data.frame(do.call(cbind, residuals_list))
  colnames(residuals_df) <- paste0(outcome_vars, "_residual")
  return(residuals_df)
}



optimal_clustering_workflow <- function(data, max_k = 8) {
  
  data_clean <- as.data.frame(scale(na.omit(data)))
  
  pca_result <- prcomp(data_clean)
  cum_var <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
  n_pcs <- which(cum_var >= 0.90)[1]
  pca_scores <- pca_result$x[, 1:n_pcs]
  
  k_range <- 2:max_k
  
  metrics <- sapply(k_range, function(k) {
    km <- kmeans(pca_scores, centers = k, nstart = 25)
    sil <- mean(silhouette(km$cluster, dist(pca_scores))[, 3])
    c(k = k, silhouette = sil)
  })

  optimal_k <- k_range[which.max(metrics["silhouette", ])]
  cat("基于轮廓系数的最佳k值:", optimal_k, "\n")

  final_clusters <- kmeans(pca_scores, centers = optimal_k, nstart = 25)

  cluster_plot <- fviz_cluster(final_clusters, data = pca_scores,
                               palette = "jco", ggtheme = theme_minimal())
  
  return(list(
    clusters = final_clusters$cluster,
    optimal_k = optimal_k,
    plot = cluster_plot,
    pca_scores = pca_scores
  ))
}


find_optimal_cuts <- function(hc, max_k = 15) {

  sil_scores <- sapply(2:max_k, function(k) {
    clusters <- cutree(hc, k)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  })

  optimal_ks <- which(diff(sign(diff(sil_scores))) == -2) + 1

  good_ks <- optimal_ks[sil_scores[optimal_ks] > 0.4]
  if (length(good_ks) == 0) {
    good_ks <- which.max(sil_scores) + 1
  }
  return(good_ks)
}


trl_plot <- function(label_text,lab_tar){
  
  library(ggplot2)
  library(dplyr)
  library(purrr)

  triangle_height <- 3.5
  triangle_base <- 6
  rect_width <- c(6.2,4,4.3)
  rect_height <- 1
  
  rect_centers <- data.frame(
    label = lab_tar,
    x = c(0, -triangle_base/2-1, triangle_base/2+1),
    y = c(triangle_height, 0, 0),
    width = rect_width,
    height = rect_height
  ) %>% mutate(
    xmin = x - width/2,
    xmax = x + width/2,
    ymin = y - height/2,
    ymax = y + height/2
  )
  rect_centers[,c(8,9)] <-  rect_centers[,c(8,9)] - 0.2
  rect_centers[,c(3)] <-  rect_centers[,c(3)] - 0.2
  calculate_edge_point <- function(from_label, to_label, rects) {
    from_rect <- filter(rects, label == from_label)
    to_rect <- filter(rects, label == to_label)
    
    # 计算方向向量
    dx <- to_rect$x - from_rect$x
    dy <- to_rect$y - from_rect$y
    len <- sqrt(dx^2 + dy^2)
    dx_unit <- dx / len
    dy_unit <- dy / len
    
    # 计算源矩形上的交点
    t_values <- c()
    if (dx_unit != 0) {
      t_values <- c(t_values, (from_rect$xmin - from_rect$x)/dx_unit, 
                    (from_rect$xmax - from_rect$x)/dx_unit)
    }
    if (dy_unit != 0) {
      t_values <- c(t_values, (from_rect$ymin - from_rect$y)/dy_unit, 
                    (from_rect$ymax - from_rect$y)/dy_unit)
    }
    
    t_min_from <- min(t_values[t_values > 0], na.rm = TRUE)
    
    t_values <- c()
    if (dx_unit != 0) {
      t_values <- c(t_values, (to_rect$xmin - to_rect$x)/(-dx_unit), 
                    (to_rect$xmax - to_rect$x)/(-dx_unit))
    }
    if (dy_unit != 0) {
      t_values <- c(t_values, (to_rect$ymin - to_rect$y)/(-dy_unit), 
                    (to_rect$ymax - to_rect$y)/(-dy_unit))
    }
    
    t_min_to <- min(t_values[t_values > 0], na.rm = TRUE)
    
    data.frame(
      from = from_label,
      to = to_label,
      x_start = from_rect$x + dx_unit * t_min_from,
      y_start = from_rect$y + dy_unit * t_min_from,
      x_end = to_rect$x - dx_unit * t_min_to,
      y_end = to_rect$y - dy_unit * t_min_to
    )
  }
  
  connections <- list(
    data.frame(from = lab_tar[2], to = lab_tar[1], color1 = "#E41A1C", color2 = "#377EB8", 
               label_text = "LKS1"),
    data.frame(from = lab_tar[1], to = lab_tar[3], color1 = "#4DAF4A", color2 = "#984EA3", 
               label_text = "LKS2"),
    data.frame(from = lab_tar[2], to = lab_tar[3], color1 = "#FF7F00", color2 = "#A65628", 
               label_text = "LKS3")
  )
  
  conn_lines <- map_df(connections, function(conn) {
    edge_points <- calculate_edge_point(from_label = conn$from, to_label = conn$to, 
                                        rects = rect_centers)
    cbind(conn, edge_points[,-c(1:2)])
  })
  
  create_arrow_lines <- function(line, offset) {
    dx <- line$x_end - line$x_start
    dy <- line$y_end - line$y_start
    len <- sqrt(dx^2 + dy^2)
    
    nx <- -dy / len * offset
    ny <- dx / len * offset
    
    arrow_obj <- arrow(angle = 20, length = unit(0.15, "inches"), type = "closed")
    
    data.frame(
      x = line$x_start + nx,
      y = line$y_start + ny,
      xend = line$x_end + nx,
      yend = line$y_end + ny,
      color = ifelse(offset > 0, line$color1, line$color2),
      group = paste0(line$from, "_", line$to, "_", sign(offset)))
  }
  
  parallel_lines <- map_df(1:nrow(conn_lines), function(i) {
    rbind(
      create_arrow_lines(conn_lines[i, ], 0.08),
      create_arrow_lines(conn_lines[i, ], -0.08)
    )
  })
  parallel_lines[c(5,6),"color"] <- parallel_lines[c(1,2),"color"]
  
  line_labels <- map_df(1:nrow(conn_lines), function(i) {
    line <- conn_lines[i, ]

    dx <- line$x_end - line$x_start
    dy <- line$y_end - line$y_start
    len <- sqrt(dx^2 + dy^2)
    
    mid_x <- (line$x_start + line$x_end) / 2
    mid_y <- (line$y_start + line$y_end) / 2
    
    nx <- -dy / len * 0.35
    ny <- dx / len * 0.35
    
    angle <- atan2(dy, dx) * 180 / pi
    if (angle > 90) angle <- angle - 180
    if (angle < -90) angle <- angle + 180
    
    data.frame(
      x = c(mid_x + nx, mid_x - nx),
      y = c(mid_y + ny, mid_y - ny),
      label = paste0(line$label_text, c("", "")),
      angle = angle,
      side = c("上", "下"),
      line_group = paste0(line$from, "_", line$to)
    )
  })
  
  line_labels <- line_labels[-4,]
  line_labels[,"label"] <- label_text 
  parallel_lines <- parallel_lines[-3,]
  line_labels[3,"y"] <-  line_labels[3,"y"] - 0.3

  fig_tem <- ggplot() +
    geom_rect(
      data = rect_centers,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
      color = "black", alpha = 0.9, show.legend = FALSE) +
    scale_fill_manual(values = c("#FFE4E1", "#E0FFFF","#F0FFF0")) +
    
    geom_text(
      data = rect_centers,
      aes(x = x, y = y, label = label),
      size = 4, fontface = "bold", color = "darkblue") + 
    geom_segment(
      data = parallel_lines,
      aes(x = x, y = y, xend = xend, yend = yend, color = color),
      linewidth = 0.8, lineend = "round", 
      arrow = arrow(angle = 10, length = unit(0.1, "inches"), type = "closed")) +
    scale_color_identity() +  
    
    geom_text(
      data = line_labels,
      aes(x = x, y = y, label = label, angle = angle),parse = TRUE,
      size = 4, fontface = "italic", color = parallel_lines$color, 
      vjust = 0.5, hjust = 0.5) +
    #annotate("text", x = -triangle_base/2, y = -1, label = "L2→L1", color = "#E41A1C", size = 4) +
    #annotate("text", x = triangle_base/2, y = -1, label = "L1→L3", color = "#4DAF4A", size = 4) +
    #annotate("text", x = 0, y = -1.5, label = "L2→L3", color = "#FF7F00", size = 4) +
    
    coord_equal(xlim = c(-triangle_base-0.5, triangle_base-0.3), 
                ylim = c(-0.5, triangle_height + 0.5)) +
    labs(title = "") +
    theme_void() +
    theme(
      plot.margin = margin(-2, -2, -1, -2),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )
  return(fig_tem)
}



enhanced_classification <- function(protein_name) {
  case_when(
    stri_detect_fixed(tolower(protein_name), "receptor") ~ "Receptor",
    stri_detect_fixed(tolower(protein_name), "kinase") ~ "Kinase",
    stri_detect_fixed(tolower(protein_name), "growth factor|differentiation") ~ "Growth_Factor",
    stri_detect_fixed(tolower(protein_name), "enzyme|dehydrogenase") ~ "Enzyme",
    stri_detect_fixed(tolower(protein_name), "adhesion|cadherin") ~ "Cell_Adhesion",
    stri_detect_fixed(tolower(protein_name), "chemokine|interleukin|tumor necrosis") ~ "Cytokine_Chemokine",
    stri_detect_fixed(tolower(protein_name), "binding protein") ~ "Binding_Protein",
    stri_detect_fixed(tolower(protein_name), "transporter") ~ "Transporter",
    TRUE ~ "Other"
  )
}
# Calculate tissue specificity function
calculate_specificity <- function(data) {
  data %>%
    group_by(name_full) %>%
    mutate(
      max_imp = max(Mean_Importance),
      specificity = ifelse(n_distinct(var) == 1, 1, 
                           max_imp / sum(Mean_Importance))
    ) %>%
    pull(specificity)
}
