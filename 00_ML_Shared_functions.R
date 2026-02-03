model_abbr <- c("GLM","GLM-EN","RF","SVM","XGBoost","NN","GBM","Ensemble")
model_full_name <- c("Generalized Linear Model",
                      "GLM with Elastic Net",
                      "Random Forest","Support Vector Machine",
                      "XGBoost","Neural Network","Grandient Boosting Machine","Ensemble")


# Plot LASSO regularization path
plot_lasso_path <- function(lasso_results) {
  plot(lasso_results$cv_model$glmnet.fit,
       xvar = "lambda",
       label = TRUE,
       main = "LASSO Regression Path")
  abline(v = log(lasso_results$best_lambda), lty = 2, col = "red")
  legend("bottomright", legend = "Optimal Lambda", lty = 2, col = "red")
}


plot_manu <- function(data){
  dat_tem <- data
  dat_tem <- mutate(dat_tem,Model = factor(Model,levels = c("glm","glmnet","rf",
                                                            "svm","xgb","nnet","gbm","Ensemble"),
                                           labels = model_abbr),
                    model_type = factor(model_type,levels = c(1:3),
                                        labels = c("Pure Protein Model","Protein + Basic Info Model",
                                                   "Full Model")),
                    var = factor(var,levels = tar_var,labels = c("Pancreas PDFF",
                                                                 "Liver PDFF",
                                                                 "Muscle fat infiltration",
                                                                 "Area of pericardial fat")))
  return(dat_tem)
}

custom_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Enhanced SHAP visualization function with direction
plot_shap_summary_with_direction <- function(shap_results, model_name = "average", top_num = 20) {
  
  # Extract importance data
  if ("average" %in% names(shap_results)) {
    # Combined results from multiple models
    importance_data <- head(shap_results$average, top_num)
  } else if ("importance" %in% names(shap_results)) {
    # Single model results
    importance_data <- head(shap_results$importance, top_num)
  } else {
    stop("Invalid SHAP results structure")
  }
  importance_data$Feature <- toupper(importance_data$Feature)
  # Create the plot
  p <- ggplot(importance_data, 
              aes(x = reorder(Feature, Mean_Abs_SHAP), 
                  y = Mean_Abs_SHAP,fill = Consensus_Direction))+
                 # y = if(!"Mean_Signed_SHAP" %in% names(importance_data)) Mean_Signed_SHAP else Mean_Abs_SHAP,
                 # fill = if("Consensus_Direction" %in% names(importance_data)) Consensus_Direction else Impact_Direction)) +
    geom_col(alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(
      name = "Impact Direction",
      values = c("Positive" = "#FF6B6B", "Negative" = "#4ECDC4", "Mixed" = "#FFD700"),
      labels = c("Positive" = "Positive Impact", "Negative" = "Negative Impact", "Mixed" = "Mixed Impact")
    ) +
    scale_x_discrete(expand = c(0,0))+
    coord_flip() +
    geom_errorbar(aes(ymin = pmax(Mean_Abs_SHAP - SD_Abs_SHAP, 0), 
                      ymax = Mean_Abs_SHAP + SD_Abs_SHAP), 
                  width = 0.2) +
    labs(
      #title = paste("SHAP Summary -", model_name),
      #subtitle = "Features sorted by importance, colored by impact direction",
      x = "Features",
      y = "Mean |SHAP|"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      axis.text.y = element_text(size = 10),
      #plot.margin = margin(5, 10, 5, 5, unit = "mm"),
      panel.grid.major.y = element_blank()
    )
  return(p)
}


# Function to create SHAP beeswarm plot with direction
plot_shap_beeswarm_detailed <- function(shap_results, model_name, top_n = 15) {
  
  # Get the specific model's SHAP results
  if ("by_model" %in% names(shap_results) && model_name %in% names(shap_results$by_model)) {
    model_results <- shap_results$by_model[[model_name]]
  } else {
    model_results <- shap_results
  }
  
  # Extract raw SHAP values
  shap_matrix <- model_results$raw_shap_matrix
  if (is.null(shap_matrix)) {
    stop("Raw SHAP matrix not available")
  }
  
  # Convert to long format
  shap_long <- as.data.frame(shap_matrix) %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Feature", values_to = "SHAP_Value") %>%
    mutate(Direction = ifelse(SHAP_Value > 0, "Positive", "Negative"))
  
  # Get top features by mean absolute SHAP
  top_features <- shap_long %>%
    group_by(Feature) %>%
    summarise(Mean_Abs_SHAP = mean(abs(SHAP_Value))) %>%
    arrange(desc(Mean_Abs_SHAP)) %>%
    head(top_n) %>%
    pull(Feature)
  
  # Filter data for top features
  plot_data <- shap_long %>%
    filter(Feature %in% top_features) %>%
    mutate(Feature = factor(Feature, levels = rev(top_features)))
  
  # Create beeswarm plot
  p <- ggplot(plot_data, aes(x = SHAP_Value, y = Feature)) +
    geom_jitter(aes(color = Direction), alpha = 0.6, size = 1, height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 18) +
    scale_color_manual(
      name = "Impact Direction",
      values = c("Positive" = "#FF6B6B", "Negative" = "#4ECDC4")
    ) +
    labs(
      title = paste("SHAP Beeswarm Plot -", model_name),
      subtitle = "Each point represents one sample, diamond shows mean SHAP value",
      x = "SHAP Value (Impact on Prediction)",
      y = "Features"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}

# Function to compare SHAP directions across models
plot_shap_direction_consistency <- function(shap_results, top_n = 20) {
  
  if (!"average" %in% names(shap_results)) {
    stop("This function requires combined SHAP results from multiple models")
  }
  
  plot_data <- head(shap_results$average, top_n)
  
  p <- ggplot(plot_data, 
              aes(x = reorder(Feature, Mean_Abs_SHAP), 
                  y = Mean_Signed_SHAP,
                  fill = Consensus_Direction,
                  alpha = Direction_Consistency)) +
    geom_col() +
    geom_errorbar(aes(ymin = Mean_Signed_SHAP - SD_Signed_SHAP, 
                      ymax = Mean_Signed_SHAP + SD_Signed_SHAP), 
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(
      name = "Consensus Direction",
      values = c("Positive" = "#FF6B6B", "Negative" = "#4ECDC4", "Mixed" = "#FFD700")
    ) +
    scale_alpha_continuous(
      name = "Direction Consistency",
      range = c(0.6, 1)
    ) +
    coord_flip() +
    labs(
      title = "SHAP Direction Consistency Across Models",
      subtitle = "Bar height shows mean SHAP value, error bars show standard deviation",
      x = "Features",
      y = "Mean Signed SHAP Value"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.major.y = element_blank()
    )
  return(p)
}


# Create SHAP summary plot
plot_shap_summary <- function(shap_results, model_name, top_n = 20) {
  
  if (!model_name %in% names(shap_results$by_model)) {
    stop("Model not found in SHAP results")
  }
  
  model_shap <- shap_results$by_model[[model_name]]
  top_features <- head(model_shap$importance, top_n)
  
  p <- ggplot(top_features, aes(x = reorder(Feature, Mean_SHAP), y = Mean_SHAP)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(
      title = paste("SHAP Feature Importance -", model_name),
      x = "Features",
      y = "Mean |SHAP value|"
    ) + theme_minimal()
  
  return(p)
}


preprocess_data <- function(data, 
                            target_var = NULL,
                            method_numeric = c("center", "scale"),
                            method_categorical = "dummy", 
                            remove_nzv = TRUE,
                            nzv_freqCut = 95/5,
                            nzv_uniqueCut = 10,
                            handle_missing = "median",  # New: missing value handling
                            correlation_threshold = 0.95) {  # New: high correlation handling
  
  # Identify target variable
  if (is.null(target_var)) {
    target_var <- names(data)[1]
    features <- data[, -1, drop = FALSE]
  } else {
    features <- data[, setdiff(names(data), target_var), drop = FALSE]
  }
  
  # New: Missing value handling
  if (handle_missing == "median") {
    features <- features %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))
  } else if (handle_missing == "knn") {
    # Use KNN imputation
    features <- VIM::kNN(features)$imp
  } else if (handle_missing == "mean") {
    features <- features %>%
      mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  } else if (handle_missing == "mode") {
    # For categorical variables
    features <- features %>%
      mutate(across(where(is.factor), ~ifelse(is.na(.), names(sort(table(.), decreasing = TRUE))[1], .)))
  }
  
  # Separate variable types
  numeric_vars <- names(features)[sapply(features, is.numeric)]
  categorical_vars <- names(features)[sapply(features, function(x) 
    is.factor(x) | is.character(x) | is.logical(x))]
  
  processed_parts <- list()
  
  # Process numeric variables
  if (length(numeric_vars) > 0) {
    numeric_data <- features[, numeric_vars, drop = FALSE]
    
    if (remove_nzv) {
      nzv_indices <- nearZeroVar(numeric_data, freqCut = nzv_freqCut, uniqueCut = nzv_uniqueCut)
      if (length(nzv_indices) > 0) {
        numeric_data <- numeric_data[, -nzv_indices, drop = FALSE]
        cat("Removed", length(nzv_indices), "near-zero variance numeric variables\n")
      }
    }
    
    if (ncol(numeric_data) > 0) {
      preprocess_params <- preProcess(numeric_data, method = method_numeric)
      processed_parts$numeric <- predict(preprocess_params, numeric_data)
    }
  }
  
  # Process categorical variables
  if (length(categorical_vars) > 0) {
    categorical_data <- features[, categorical_vars, drop = FALSE]
    
    # Convert to factors
    categorical_data <- categorical_data %>%
      mutate(across(where(is.character), as.factor))
    
    if (remove_nzv) {
      nzv_indices <- nearZeroVar(categorical_data, freqCut = nzv_freqCut, uniqueCut = nzv_uniqueCut)
      if (length(nzv_indices) > 0) {
        categorical_data <- categorical_data[, -nzv_indices, drop = FALSE]
        cat("Removed", length(nzv_indices), "near-zero variance categorical variables\n")
      }
    }
    
    if (ncol(categorical_data) > 0) {
      if (method_categorical == "dummy") {
        # One-hot encoding
        dummies <- dummyVars(~ ., data = categorical_data, fullRank = TRUE)
        processed_parts$categorical <- predict(dummies, categorical_data) %>% as.data.frame()
      } else if (method_categorical == "label") {
        # Label encoding (for tree-based models)
        processed_parts$categorical <- categorical_data %>%
          mutate(across(everything(), as.numeric))
      } else {
        stop("Unknown categorical processing method. Use 'dummy' or 'label'")
      }
    }
  }
  
  # Combine all processed parts
  if (length(processed_parts) == 2) {
    processed_features <- cbind(processed_parts$numeric, processed_parts$categorical)
  } else if (length(processed_parts) == 1) {
    processed_features <- processed_parts[[1]]
  } else {
    stop("No features remaining after preprocessing")
  }
  
  # New: High correlation feature handling
  if (correlation_threshold < 1) {
    numeric_features <- processed_features[, sapply(processed_features, is.numeric), drop = FALSE]
    if (ncol(numeric_features) > 1) {
      cor_matrix <- cor(numeric_features, use = "complete.obs")
      high_cor <- findCorrelation(cor_matrix, cutoff = correlation_threshold)
      if (length(high_cor) > 0) {
        processed_features <- processed_features[, -high_cor, drop = FALSE]
        cat("Removed", length(high_cor), "highly correlated features\n")
      }
    }
  }
  
  # Final data frame
  result <- cbind(data[, target_var, drop = FALSE], processed_features)
  names(result)[1] <- target_var
  return(result)
}


analyze_nonlinear_effects <- function(selected_features) {
  nonlinear_effects <- list()
  
  for(feature in selected_features) {
    if(grepl("_spline", feature)) {
      
      base_var <- gsub("_spline[0-9]+", "", feature)
      if(!base_var %in% names(nonlinear_effects)) {
        nonlinear_effects[[base_var]] <- c()
      }
      nonlinear_effects[[base_var]] <- c(nonlinear_effects[[base_var]], feature)
    }
  }
  strong_nonlinear <- names(nonlinear_effects)[sapply(nonlinear_effects, length) > 1]
  
  return(list(
    all_nonlinear = names(nonlinear_effects),
    strong_nonlinear = strong_nonlinear
  ))
}


create_direction_heatmap <- function(model_level_shap) {
  
  model_level_shap <- mutate(model_level_shap, 
                             Feature = toupper(Feature),
                             Model = factor(Model,levels = c("glm","glmnet","rf",
                                                             "svm","xgb","nnet","gbm"),
                                            labels = model_abbr[-8]))

  consistency_df <- model_level_shap %>%
    group_by(Feature) %>%
    summarise(
      total_models = n(),
      positive_count = sum(Impact_Direction == "Positive", na.rm = TRUE),
      negative_count = sum(Impact_Direction == "Negative", na.rm = TRUE),
      consistency = max(positive_count, negative_count) / total_models,
      dominant_direction = ifelse(positive_count >= negative_count, "Positive", "Negative")
    ) %>%
    arrange(desc(consistency), dominant_direction) %>%
    mutate(Feature_ordered = factor(Feature, levels = unique(Feature)))
  
  model_level_shap <- model_level_shap %>%
    mutate(Feature = factor(Feature, levels = rev(levels(consistency_df$Feature_ordered))))
  
  heatmap_plot <- ggplot(model_level_shap, 
                         aes(x = Model, y = Feature, fill = Impact_Direction)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(
      values = c("Positive" = "#FF9AA2", "Negative" = "#89CFF0"),
      na.value = "gray90") +
    labs(
      x = "Machine Learning Model",
      y = "Protein Features",
      #title = "SHAP Value Directions Across Models",
      #subtitle = "Features sorted by directional consistency across models",
      fill = "Direction") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      plot.subtitle = element_text(size = 10, color = "gray40"))
  
  return(heatmap_plot)
}



extract_all_importance <- function(shap_results) {
  importance_list <- list()
  
  for (model_name in names(shap_results$by_model)) {
    importance_df <- shap_results$by_model[[model_name]]$importance
    importance_df$Model <- model_name
    importance_list[[model_name]] <- importance_df
  }
  
  combined_importance <- do.call(rbind, importance_list)
  rownames(combined_importance) <- NULL
  
  return(list(
    by_model = importance_list,
    combined = combined_importance
  ))
}

# Analyze variable importance across models
analyze_variable_importance <- function(models, top_n = 15) {
  importance_list <- list()
  
  # Extract variable importance from each model
  for (model_name in names(models)) {
    imp <- tryCatch({
      varImp(models[[model_name]])$importance
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(imp)) {
      if ("Overall" %in% colnames(imp)) {
        imp_df <- data.frame(
          Protein = rownames(imp),
          Importance = imp$Overall,
          Model = model_name
        )
      } else {
        # For multi-column importance, take first column or mean
        if (ncol(imp) == 1) {
          imp_df <- data.frame(
            Protein = rownames(imp),
            Importance = imp[, 1],
            Model = model_name
          )
        } else {
          imp_mean <- rowMeans(imp, na.rm = TRUE)
          imp_df <- data.frame(
            Protein = rownames(imp),
            Importance = imp_mean,
            Model = model_name
          )
        }
      }
      importance_list[[model_name]] <- imp_df
    }
  }
  
  # Combine all importance scores
  all_imp <- bind_rows(importance_list)
  
  # Calculate average importance across models
  avg_importance <- all_imp %>%
    group_by(Protein) %>%
    summarise(
      Mean_Importance = mean(Importance, na.rm = TRUE),
      SD_Importance = sd(Importance, na.rm = TRUE),
      N_Models = n()
    ) %>%
    arrange(desc(Mean_Importance)) %>%
    head(top_n)
  
  return(list(
    detailed = all_imp,
    average = avg_importance
  ))
}
