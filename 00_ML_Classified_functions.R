
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

# Function to analyze gender distribution of binary outcomes
analyze_gender_distribution <- function(data, binary_vars) {
  
  results <- list()
  
  for (var in binary_vars) {
    # Create contingency table
    cont_table <- table(data$sex, data[[var]])
    
    # Calculate proportions
    prop_table <- prop.table(cont_table, margin = 1)
    
    # Statistical test (Chi-square)
    chi_test <- chisq.test(cont_table)
    
    results[[var]] <- list(
      contingency_table = cont_table,
      proportion_table = prop_table,
      chi_square_test = chi_test
    )
    
    cat("=== Analysis for", var, "===\n")
    cat("Contingency Table:\n")
    print(cont_table)
    cat("\nProportions by Gender:\n")
    print(prop_table)
    cat("\nChi-square test p-value:", chi_test$p.value, "\n\n")
  }
  
  return(results)
}

# Function to validate gender-specific thresholds
validate_gender_thresholds <- function(data, thresholds) {
  
  validation_results <- list()
  
  for (gender in c("Male", "Female")) {
    gender_code <- ifelse(gender == "Male", 1, 0)
    gender_data <- data %>% filter(sex == gender_code)
    
    cat("=== Validation for", gender, "===\n")
    
    for (measurement in names(thresholds[[gender]])) {
      threshold <- thresholds[[gender]][[measurement]]
      prevalence <- mean(gender_data[[measurement]] > threshold, na.rm = TRUE)
      
      cat(sprintf("%s: Threshold = %.1f, Prevalence = %.1f%%\n", 
                  measurement, threshold, prevalence * 100))
      
      validation_results[[paste0(gender, "_", measurement)]] <- list(
        threshold = threshold,
        prevalence = prevalence,
        n_above_threshold = sum(gender_data[[measurement]] > threshold, na.rm = TRUE),
        n_total = sum(!is.na(gender_data[[measurement]]))
      )
    }
    cat("\n")
  }
  
  return(validation_results)
}

# Function to find optimal gender-specific thresholds (if needed)
find_optimal_gender_thresholds <- function(data, method = "youden") {
  
  library(pROC)
  
  optimal_thresholds <- list(Male = list(), Female = list())
  
  for (gender in c("Male", "Female")) {
    gender_code <- ifelse(gender == "Male", 1, 0)
    gender_data <- data %>% filter(sex == gender_code)
    
    cat("=== Finding optimal thresholds for", gender, "===\n")
    
    for (measurement in c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction", 
                          "Muscle_fat_infiltration", "area_of_pericardial_fat")) {
      
      # This is a placeholder - in real analysis, you would have a gold standard
      # For demonstration, we'll use median as a simple cutoff
      current_data <- gender_data[[measurement]]
      optimal_threshold <- median(current_data, na.rm = TRUE)
      
      optimal_thresholds[[gender]][[measurement]] <- optimal_threshold
      
      cat(sprintf("%s: Optimal threshold = %.2f\n", measurement, optimal_threshold))
    }
    cat("\n")
  }
  
  return(optimal_thresholds)
}



# New helper functions for binary classification
compare_three_classification_models <- function(models_results) {
  
  comparison_df <- data.frame()
  
  for (i in 1:length(models_results)) {
    model <- models_results[[i]]
    best_performance <- model$model_performance[which.max(model$model_performance$AUC), ]
    
    comparison_df <- rbind(comparison_df, data.frame(
      Model_Type = model$model_type,
      Model_Description = model$model_description,
      Best_Model = rownames(best_performance),
      AUC = best_performance$AUC,
      Accuracy = best_performance$Accuracy,
      Sensitivity = best_performance$Sensitivity,
      Specificity = best_performance$Specificity,
      Feature_Count = model$feature_count,
      Protein_Score_AUC_train = if(!is.null(model$protein_score_performance$train)) model$protein_score_performance$train$auc else NA,
      Protein_Score_AUC_test = if(!is.null(model$protein_score_performance$test)) model$protein_score_performance$test$auc else NA
    ))
  }
  return(comparison_df)
}

plot_classification_model_comparison <- function(comparison_df) {
  # Create performance comparison plot
  p1 <- ggplot(comparison_df, aes(x = Model_Description, y = AUC, fill = Model_Description)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5) +
    labs(title = "Model Performance Comparison (AUC)",
         x = "Model Type", y = "AUC") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(comparison_df, aes(x = Model_Description, y = Feature_Count, fill = Model_Description)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Feature_Count), vjust = -0.5) +
    labs(title = "Number of Selected Features",
         x = "Model Type", y = "Feature Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(auc_plot = p1, feature_plot = p2))
}



print_classification_model_comparison_summary <- function(models_results, target_var) {
  cat("\n=== Classification Model Comparison Summary - Target Variable:", target_var, "===\n")
  
  for (i in 1:length(models_results)) {
    model <- models_results[[i]]
    best_perf <- model$model_performance[which.max(model$model_performance$AUC), ]
    
    cat(sprintf("\nModel %d (%s):\n", i, model$model_description))
    cat(sprintf("  Best Model: %s\n", rownames(best_perf)))
    cat(sprintf("  AUC: %.4f\n", best_perf$AUC))
    cat(sprintf("  Accuracy: %.4f\n", best_perf$Accuracy))
    cat(sprintf("  Sensitivity: %.4f\n", best_perf$Sensitivity))
    cat(sprintf("  Specificity: %.4f\n", best_perf$Specificity))
    cat(sprintf("  Feature Count: %d\n", model$feature_count))
    
    if(!is.null(model$protein_score_performance$test)) {
      cat(sprintf("  Protein Score AUC (Test Set): %.4f\n", model$protein_score_performance$test$auc))
    }
    
    # Display top 5 most important features
    if(!is.null(model$variable_importance$average)) {
      top_features <- head(model$variable_importance$average, 5)
      cat("  Top 5 Most Important Features:\n")
      for (j in 1:nrow(top_features)) {
        cat(sprintf("    %d. %s (Importance: %.3f)\n", 
                    j, top_features$Feature[j], top_features$Mean_Importance[j]))
      }
    }
  }
  
  # Performance improvement analysis
  if (length(models_results) >= 2) {
    auc_improvement_1_to_2 <- models_results[[2]]$model_performance$AUC[1] - models_results[[1]]$model_performance$AUC[1]
    auc_improvement_2_to_3 <- models_results[[3]]$model_performance$AUC[1] - models_results[[2]]$model_performance$AUC[1]
    
    cat(sprintf("\nPerformance Improvement Analysis:\n"))
    cat(sprintf("  AUC Improvement from Adding Basic Info: %.4f\n", auc_improvement_1_to_2))
    cat(sprintf("  AUC Improvement from Adding Body Composition: %.4f\n", auc_improvement_2_to_3))
    cat(sprintf("  Total AUC Improvement: %.4f\n", auc_improvement_1_to_2 + auc_improvement_2_to_3))
  }
}

# LASSO feature selection for classification
perform_lasso_selection_classification <- function(train_data, alpha = 1) {
  
  x <- as.matrix(train_data[, -1])
  y <- train_data$target
  
  # Ensure y is numeric for glmnet (0/1)
  if(is.factor(y)) {
    y_numeric <- as.numeric(y) - 1
  } else {
    y_numeric <- y
  }
  
  # Cross-validated LASSO for classification
  cv_lasso <- cv.glmnet(x, y_numeric, 
                        family = "binomial",  # Changed to binomial for classification
                        alpha = alpha,
                        nfolds = 10,
                        type.measure = "auc")  # Use AUC for classification
  
  best_lambda <- cv_lasso$lambda.1se
  
  # Fit final LASSO model
  lasso_model <- glmnet(x, y_numeric, 
                        family = "binomial",  # Changed to binomial
                        alpha = alpha, 
                        lambda = best_lambda)
  
  # Extract coefficients and selected features
  lasso_coef <- as.matrix(coef(lasso_model))
  selected_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
  selected_features <- selected_features[selected_features != "(Intercept)"]
  
  coef_values <- lasso_coef[which(lasso_coef != 0), ]
  names(coef_values) <- rownames(lasso_coef)[which(lasso_coef != 0)]
  coef_values <- coef_values[names(coef_values) != "(Intercept)"]
  
  results <- list(
    cv_model = cv_lasso,
    final_model = lasso_model,
    selected_features = selected_features,
    coefficients = coef_values,
    best_lambda = best_lambda
  )
  return(results)
}

# Compare classification model performance
compare_classification_models <- function(models, test_data) {
  performance_df <- data.frame()
  predictions_list <- list()
  
  for (model_name in names(models)) {
    # Generate predictions
    pred_prob <- predict(models[[model_name]], test_data, type = "prob")
    pred_class <- predict(models[[model_name]], test_data)
    true_values <- test_data$target
    
    # Calculate performance metrics
    cm <- confusionMatrix(pred_class, true_values)
    auc_val <- roc(true_values, pred_prob[,2])$auc
    
    performance_df <- rbind(performance_df, data.frame(
      Model = model_name,
      AUC = as.numeric(auc_val),
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      Kappa = cm$overall["Kappa"]
    ))
    
    predictions_list[[model_name]] <- list(
      probabilities = pred_prob,
      classes = pred_class
    )
  }
  
  performance_df <- performance_df[order(-performance_df$AUC), ]
  return(list(
    performance = performance_df,
    predictions = predictions_list
  ))
}

# Create ensemble for classification
create_classification_ensemble <- function(predictions_list, true_values, method = "weighted") {
  
  # Extract probabilities for positive class
  prob_list <- lapply(predictions_list, function(x) x$probabilities[,2])
  prob_matrix <- do.call(cbind, prob_list)
  
  if (method == "weighted") {
    # Weighted average based on AUC performance
    weights <- sapply(predictions_list, function(x) {
      roc_obj <- roc(true_values, x$probabilities[,2])
      auc(roc_obj)
    })
    weights <- weights / sum(weights)
    
    ensemble_prob <- rowSums(prob_matrix * weights)
  } else {
    # Simple average
    ensemble_prob <- rowMeans(prob_matrix)
  }
  
  # Convert to class predictions
  positive_class <- levels(true_values)[2]
  negative_class <- levels(true_values)[1]
  ensemble_class <- factor(ifelse(ensemble_prob > 0.5, positive_class, negative_class),
                           levels = levels(true_values))
  
  # Calculate ensemble performance
  cm_ensemble <- confusionMatrix(ensemble_class, true_values)
  auc_ensemble <- roc(true_values, ensemble_prob)$auc
  
  ensemble_performance <- data.frame(
    Model = "Ensemble",
    AUC = as.numeric(auc_ensemble),
    Accuracy = cm_ensemble$overall["Accuracy"],
    Sensitivity = cm_ensemble$byClass["Sensitivity"],
    Specificity = cm_ensemble$byClass["Specificity"],
    Kappa = cm_ensemble$overall["Kappa"]
  )
  
  return(list(
    performance = ensemble_performance,
    predictions = list(
      probabilities = data.frame(
        Low = 1 - ensemble_prob, 
        High = ensemble_prob
      ),
      classes = ensemble_class
    )
  ))
}

# Build protein score for classification
build_protein_score_classification <- function(data, coefficients) {
  
  # Match available proteins with coefficient names
  available_proteins <- names(coefficients)[names(coefficients) %in% colnames(data)]
  
  if (length(available_proteins) == 0) {
    stop("No matching proteins found for score calculation")
  }
  
  selected_coef <- coefficients[available_proteins]
  selected_data <- data[, available_proteins, drop = FALSE]
  
  # Calculate protein score: linear combination of protein levels
  protein_score <- as.matrix(selected_data) %*% selected_coef
  
  return(as.numeric(protein_score))
}

# Evaluate protein score performance for classification
evaluate_protein_score_classification <- function(protein_score, true_values, 
                                                  dataset_name = "") {
  
  # Calculate ROC and AUC
  roc_obj <- roc(true_values, protein_score)
  auc_val <- pROC::auc(roc_obj)
  
  # Find optimal cutoff
  coords_obj <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  optimal_cutoff <- coords_obj["threshold"]
  
  # Calculate predictions at optimal cutoff
  pred_class <- ifelse(protein_score > as.numeric(optimal_cutoff), levels(true_values)[2], levels(true_values)[1])
  pred_class <- factor(pred_class, levels = levels(true_values))
  
  # Calculate confusion matrix metrics
  cm <- confusionMatrix(pred_class, true_values)
  
  results <- list(
    auc = auc_val,
    optimal_cutoff = optimal_cutoff,
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"],
    accuracy = cm$overall["Accuracy"],
    roc_obj = roc_obj,
    plot_data = data.frame(
      True_Value = true_values,
      Protein_Score = protein_score
    )
  )
  return(results)
}

# Create classification diagnostic plots
create_classification_plots <- function(train_results, test_results) {
  
  # ROC curve plot
  roc_data <- data.frame(
    FPR = c(1 - train_results$roc_obj$specificities, 1 - test_results$roc_obj$specificities),
    TPR = c(train_results$roc_obj$sensitivities, test_results$roc_obj$sensitivities),
    Dataset = c(rep(paste0("Training (AUC = ", round(train_results$auc, 3), ")"), 
                    length(train_results$roc_obj$sensitivities)),
                rep(paste0("Testing (AUC = ", round(test_results$auc, 3), ")"), 
                    length(test_results$roc_obj$sensitivities)))
  )
  
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Dataset)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = "ROC Curves for Protein Score",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Score distribution plot
  train_plot_data <- train_results$plot_data
  train_plot_data$Dataset <- "Training"
  test_plot_data <- test_results$plot_data
  test_plot_data$Dataset <- "Testing"
  all_plot_data <- rbind(train_plot_data, test_plot_data)
  
  dist_plot <- ggplot(all_plot_data, aes(x = Protein_Score, fill = True_Value)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ Dataset, ncol = 2) +
    labs(title = "Protein Score Distribution by True Class",
         x = "Protein Score",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(
    roc_plot = roc_plot,
    distribution_plot = dist_plot
  ))
}

# LASSO for fixed features in classification
lasso_fixed_features_classification <- function(data, features, alpha = 1, nfolds = 10) {
  
  # Ensure only pre-selected features are used
  X_subset <- data[, features]
  
  # Convert target to numeric if factor
  y <- data$target
  if(is.factor(y)) {
    y_numeric <- as.numeric(y) - 1
  } else {
    y_numeric <- y
  }
  
  # Use cross-validation to select lambda
  cv_lasso <- cv.glmnet(
    x = as.matrix(X_subset),
    y = y_numeric,
    family = "binomial",  # Changed to binomial
    alpha = alpha,
    nfolds = nfolds,
    type.measure = "auc",  # Use AUC for classification
    lambda.min.ratio = 0.001
  )
  
  # Use lambda.min to retain more coefficients
  coef_matrix <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  coefficients <- coef_matrix[-1, 1]  # Remove intercept
  names(coefficients) <- features
  
  return(coefficients)
}

# Find the best classification model
find_best_classification_model <- function(models_list) {
  best_auc <- -Inf
  best_model_index <- 1
  best_model_name <- ""
  
  for (i in 1:length(models_list)) {
    model_perf <- models_list[[i]]$model_performance
    current_best_auc <- max(model_perf$AUC, na.rm = TRUE)
    
    if (current_best_auc > best_auc) {
      best_auc = current_best_auc
      best_model_index = i
      best_model_name = model_perf$Model[which.max(model_perf$AUC)]
    }
  }
  
  return(list(
    best_model_index = best_model_index,
    best_model_name = best_model_name,
    best_auc = best_auc
  ))
}

# Print hierarchical classification model summary
print_hierarchical_classification_summary <- function(models_list, target_var) {
  cat("\n")
  cat("===========================================\n")
  cat("Classification Target Variable:", target_var, "\n")
  cat("===========================================\n")
  
  # Protein selection summary
  cat("Protein Selection Results:\n")
  cat("- Number of selected proteins:", length(models_list[[1]]$selected_proteins), "\n")
  cat("- Selected proteins:", paste(models_list[[1]]$selected_proteins, collapse = ", "), "\n")
  
  # Model performance summary
  cat("\nModel Performance Comparison:\n")
  for (i in 1:length(models_list)) {
    best_perf <- models_list[[i]]$model_performance[which.max(models_list[[i]]$model_performance$AUC), ]
    cat(sprintf("- %s: AUC = %.3f, Accuracy = %.3f, Feature Count = %d\n",
                models_list[[i]]$model_description,
                best_perf$AUC,
                best_perf$Accuracy,
                models_list[[i]]$feature_count))
  }
  
  cat("===========================================\n\n")
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


# Prepare LASSO data for classification
prepare_lasso_data_classification <- function(data, predictors, spline_df = 3) {
  X_raw <- data[, predictors, drop = FALSE]
  
  spline_features <- data.frame(row.names = rownames(data))
  for(pred in predictors) {
    if(is.numeric(X_raw[[pred]])) {
      spline_basis <- ns(X_raw[[pred]], df = spline_df)
      colnames(spline_basis) <- paste0(pred, "_spline", 1:spline_df)
      spline_features <- cbind(spline_features, spline_basis)
    }
  }
  X_final <- as.matrix(spline_features)
  return(list(X = X_final, feature_names = colnames(X_final)))
}



# Train multiple classification models
train_classification_models <- function(train_data) {
  
  models <- list()
  
  # Ensure target is factor with proper levels
  if(!is.factor(train_data$target)) {
    train_data$target <- as.factor(train_data$target)
  }
  
  # Configure cross-validation for classification
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    savePredictions = "final",
    classProbs = TRUE,           # Important for classification
    summaryFunction = twoClassSummary,  # For binary classification metrics
    sampling = "up"              # Upsampling to handle class imbalance
  )
  
  # 1. Logistic Regression
  cat("Training Logistic Regression Model...\n")
  models$glm <- train(
    target ~ .,
    data = train_data,
    trControl = ctrl,
    method = "glm",
    family = binomial(link = "logit"),  # Changed to binomial for classification
    metric = "ROC"                     # Use ROC/AUC as primary metric
  )
  
  # 2. Random Forest
  cat("Training Random Forest Model...\n")
  models$rf <- train(
    target ~ .,
    data = train_data,
    method = "rf",
    trControl = ctrl,
    metric = "ROC",                    # Use ROC for classification
    tuneLength = 5,
    importance = TRUE
  )
  
  # 3. XGBoost for Classification
  cat("Training XGBoost Model...\n")
  models$xgb <- train(
    target ~ .,
    data = train_data,
    method = "xgbTree",
    trControl = ctrl,
    metric = "ROC",                    # Use ROC for classification
    tuneLength = 5
    # Note: xgboost automatically detects binary classification
  )
  
  # 4. Elastic Net for Logistic Regression
  cat("Training Elastic Net (Logistic) Model...\n")
  models$glmnet <- train(
    target ~ .,
    data = train_data,
    method = "glmnet",
    trControl = ctrl,
    metric = "ROC",                    # Use ROC for classification
    tuneLength = 5,
    family = "binomial"               # Specify binomial family for classification
  )
  
  # 5. Support Vector Machine with Radial Basis Function
  cat("Training Support Vector Machine Model...\n")
  models$svm <- train(
    target ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = ctrl,
    metric = "ROC",                    # Use ROC for classification
    tuneLength = 5,
    probability = TRUE                # Enable probability predictions
  )
  
  # 6. Neural Network for Classification
  cat("Training Neural Network...\n")
  models$nnet <- train(
    target ~ .,
    data = train_data,
    method = "nnet",
    trControl = ctrl,
    metric = "ROC",                    # Use ROC for classification
    tuneLength = 5,
    trace = FALSE,
    tuneGrid = expand.grid(
      size = c(1,3,5,7,9),
      decay = c(0.0001,0.001,0.01,0.1,0.5)
    ),
    linout = FALSE,                   # Changed to FALSE for classification
    MaxNWts = 1000,                   # Increased weight limit
    maxit = 200
  )
  
  # 7. Additional: Gradient Boosting Machine
  cat("Training Gradient Boosting Machine...\n")
  models$gbm <- train(
    target ~ .,
    data = train_data,
    method = "gbm",
    trControl = ctrl,
    metric = "ROC",
    tuneLength = 5,
    verbose = FALSE
  )
  
  return(models)
}



# Model comparison function for classification
compare_classification_models <- function(models, test_data) {
  
  performance_list <- list()
  predictions_list <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Get predictions
    pred_prob <- predict(model, test_data, type = "prob")
    pred_class <- predict(model, test_data)
    
    # Calculate performance metrics
    cm <- confusionMatrix(pred_class, test_data$target)
    auc_val <- roc(test_data$target, pred_prob[,2])$auc
    
    performance_list[[model_name]] <- data.frame(
      Model = model_name,
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      AUC = as.numeric(auc_val),
      Kappa = cm$overall["Kappa"]
    )
    
    predictions_list[[model_name]] <- list(
      probabilities = pred_prob,
      classes = pred_class
    )
  }
  
  performance_df <- do.call(rbind, performance_list)
  rownames(performance_df) <- NULL
  
  return(list(
    performance = performance_df,
    predictions = predictions_list
  ))
}

# Ensemble method for classification
create_classification_ensemble <- function(predictions_list, test_target, method = "weighted") {
  
  # Extract probabilities for positive class
  prob_list <- lapply(predictions_list, function(x) x$probabilities[,2])
  prob_matrix <- do.call(cbind, prob_list)
  
  if (method == "weighted") {
    # Simple average ensemble
    ensemble_prob <- rowMeans(prob_matrix)
  } else if (method == "stacking") {
    # Could implement stacking with meta-learner
    ensemble_prob <- rowMeans(prob_matrix)  # Placeholder
  }
  
  # Convert to class predictions
  ensemble_class <- factor(ifelse(ensemble_prob > 0.5, 
                                  levels(test_target)[2], 
                                  levels(test_target)[1]),
                           levels = levels(test_target))
  
  # Calculate ensemble performance
  cm_ensemble <- confusionMatrix(ensemble_class, test_target)
  auc_ensemble <- roc(test_target, ensemble_prob)$auc
  
  ensemble_performance <- data.frame(
    Model = "Ensemble",
    Accuracy = cm_ensemble$overall["Accuracy"],
    Sensitivity = cm_ensemble$byClass["Sensitivity"],
    Specificity = cm_ensemble$byClass["Specificity"],
    AUC = as.numeric(auc_ensemble),
    Kappa = cm_ensemble$overall["Kappa"]
  )
  
  return(list(
    performance = ensemble_performance,
    predictions = list(
      probabilities = data.frame(Low = 1 - ensemble_prob, High = ensemble_prob),
      classes = ensemble_class
    )
  ))
}


# SHAP analysis for classification models
analyze_shap_values_classification <- function(models, test_data, top_n = 50) {
  library(fastshap)
  
  shap_results <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Prepare prediction function for SHAP (returns probability of positive class)
    pfun <- function(object, newdata) {
      if (model_name == "xgb") {
        predict(object, as.matrix(newdata), type = "prob")[,2]
      } else {
        predict(object, newdata, type = "prob")[,2]
      }
    }
    
    # Calculate SHAP values
    shap_values <- fastshap::explain(
      model, 
      X = test_data[, -which(names(test_data) == "target")],
      pred_wrapper = pfun,
      nsim = 100
    )
    
    # Calculate both absolute and signed SHAP values
    mean_abs_shap <- colMeans(abs(shap_values))
    mean_signed_shap <- colMeans(shap_values)
    
    # Create importance data frame with both magnitude and direction
    shap_importance <- data.frame(
      Feature = names(mean_abs_shap),
      Mean_Abs_SHAP = as.numeric(mean_abs_shap),
      Mean_Signed_SHAP = as.numeric(mean_signed_shap),
      Impact_Direction = ifelse(mean_signed_shap > 0, "Positive", "Negative")
    )
    
    shap_importance <- shap_importance[order(shap_importance$Mean_Abs_SHAP, decreasing = TRUE), ]
    
    shap_results[[model_name]] <- list(
      shap_values = shap_values,
      importance = head(shap_importance, top_n),
      model_type = model_name,
      raw_shap_matrix = shap_values
    )
  }
  
  # Combine SHAP importance across all models
  all_features <- unique(unlist(lapply(shap_results, function(x) x$importance$Feature)))
  combined_importance <- data.frame(
    Feature = all_features,
    Mean_Abs_SHAP = 0,
    Mean_Signed_SHAP = 0,
    SD_Abs_SHAP = 0,
    SD_Signed_SHAP = 0,
    Models_Count = 0,
    Positive_Models = 0,
    Negative_Models = 0,
    Consensus_Direction = "",
    Direction_Consistency = 0
  )
  
  for (i in 1:nrow(combined_importance)) {
    feature <- combined_importance$Feature[i]
    abs_shap_values <- c()
    signed_shap_values <- c()
    directions <- c()
    
    for (model_name in names(shap_results)) {
      model_importance <- shap_results[[model_name]]$importance
      if (feature %in% model_importance$Feature) {
        abs_val <- model_importance$Mean_Abs_SHAP[model_importance$Feature == feature]
        signed_val <- model_importance$Mean_Signed_SHAP[model_importance$Feature == feature]
        direction <- model_importance$Impact_Direction[model_importance$Feature == feature]
        
        abs_shap_values <- c(abs_shap_values, abs_val)
        signed_shap_values <- c(signed_shap_values, signed_val)
        directions <- c(directions, direction)
      }
    }
    
    if (length(abs_shap_values) > 0) {
      combined_importance$Mean_Abs_SHAP[i] <- mean(abs_shap_values)
      combined_importance$Mean_Signed_SHAP[i] <- mean(signed_shap_values)
      combined_importance$SD_Abs_SHAP[i] <- sd(abs_shap_values)
      combined_importance$SD_Signed_SHAP[i] <- sd(signed_shap_values)
      combined_importance$Models_Count[i] <- length(abs_shap_values)
      combined_importance$Positive_Models[i] <- sum(directions == "Positive")
      combined_importance$Negative_Models[i] <- sum(directions == "Negative")
      
      # Determine consensus direction
      if (combined_importance$Positive_Models[i] > combined_importance$Negative_Models[i]) {
        combined_importance$Consensus_Direction[i] <- "Positive"
      } else if (combined_importance$Negative_Models[i] > combined_importance$Positive_Models[i]) {
        combined_importance$Consensus_Direction[i] <- "Negative"
      } else {
        combined_importance$Consensus_Direction[i] <- "Mixed"
      }
      
      # Calculate direction consistency
      combined_importance$Direction_Consistency[i] <- 
        max(combined_importance$Positive_Models[i], combined_importance$Negative_Models[i]) / 
        combined_importance$Models_Count[i]
    }
  }
  
  combined_importance <- combined_importance[order(combined_importance$Mean_Abs_SHAP, decreasing = TRUE), ]
  
  return(list(
    by_model = shap_results,
    average = combined_importance,
    top_features = head(combined_importance, top_n)
  ))
}


# Print classification model summary
print_classification_model_summary <- function(all_models_results, variable_name) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("CLASSIFICATION MODEL SUMMARY:", variable_name, "\n")
  cat(strrep("=", 60), "\n\n")
  
  for (i in 1:length(all_models_results)) {
    model_info <- all_models_results[[i]]
    
    cat("Model Type:", model_info$model_description, "\n")
    cat("Number of features:", model_info$feature_count, "\n")
    cat("Selected proteins:", length(model_info$selected_proteins), "\n\n")
    
    # Print performance metrics
    cat("Performance Metrics:\n")
    print(model_info$model_performance)
    cat("\n")
    
    # Variable importance (top 5)
    if (!is.null(model_info$variable_importance$average)) {
      cat("Top 5 Most Important Features:\n")
      top_features <- head(model_info$variable_importance$average, 5)
      print(top_features[, c("Protein", "Mean_Importance")])
    }
    cat("\n", strrep("-", 40), "\n\n")
  }
}
