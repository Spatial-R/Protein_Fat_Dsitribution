
# New helper functions
compare_three_models <- function(models_results) {
  
  comparison_df <- data.frame()
  
  for (i in 1:length(models_results)) {
    model <- models_results[[i]]
    best_performance <- model$model_performance[which.max(model$model_performance$R_squared), ]
    
    comparison_df <- rbind(comparison_df, data.frame(
      Model_Type = model$model_type,
      Model_Description = model$model_description,
      Best_Model = rownames(best_performance),
      R2 = best_performance$R_squared,
      RMSE = best_performance$RMSE,
      MAE = best_performance$MAE,
      Feature_Count = model$feature_count,
      Protein_Score_R2_train = model$protein_score_performance$train$r_squared,
      Protein_Score_R2_test = model$protein_score_performance$test$r_squared
    ))
  }
  return(comparison_df)
}


plot_model_comparison <- function(comparison_df) {
  # Create performance comparison plot
  p1 <- ggplot(comparison_df, aes(x = Model_Description, y = R2, fill = Model_Description)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("%.3f", R2)), vjust = -0.5) +
    labs(title = "Model Performance Comparison (R-squared)",
         x = "Model Type", y = "R-squared") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(comparison_df, aes(x = Model_Description, y = Feature_Count, fill = Model_Description)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = Feature_Count), vjust = -0.5) +
    labs(title = "Number of Selected Features",
         x = "Model Type", y = "Feature Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(r2_plot = p1, feature_plot = p2))
}


print_model_comparison_summary <- function(models_results, target_var) {
  cat("\n=== Model Comparison Summary - Target Variable:", target_var, "===\n")
  
  for (i in 1:length(models_results)) {
    model <- models_results[[i]]
    best_perf <- model$model_performance[which.max(model$model_performance$R2), ]
    
    cat(sprintf("\nModel %d (%s):\n", i, model$model_description))
    cat(sprintf("  Best Model: %s\n", rownames(best_perf)))
    cat(sprintf("  R²: %.4f\n", best_perf$R2))
    cat(sprintf("  RMSE: %.4f\n", best_perf$RMSE))
    cat(sprintf("  Feature Count: %d\n", model$feature_count))
    cat(sprintf("  Protein Score R² (Test Set): %.4f\n", model$protein_score_performance$test$r_squared))
    
    # Display top 5 most important features
    top_features <- head(model$variable_importance$average, 5)
    cat("  Top 5 Most Important Features:\n")
    for (j in 1:nrow(top_features)) {
      cat(sprintf("    %d. %s (Importance: %.3f)\n", 
                  j, top_features$Protein[j], top_features$Mean_Importance[j]))
    }
  }
  
  # Performance improvement analysis
  if (length(models_results) >= 2) {
    r2_improvement_1_to_2 <- models_results[[2]]$model_performance$R2[1] - models_results[[1]]$model_performance$R2[1]
    r2_improvement_2_to_3 <- models_results[[3]]$model_performance$R2[1] - models_results[[2]]$model_performance$R2[1]
    
    cat(sprintf("\nPerformance Improvement Analysis:\n"))
    cat(sprintf("  R² Improvement from Adding Basic Info: %.4f\n", r2_improvement_1_to_2))
    cat(sprintf("  R² Improvement from Adding Body Composition: %.4f\n", r2_improvement_2_to_3))
    cat(sprintf("  Total R² Improvement: %.4f\n", r2_improvement_1_to_2 + r2_improvement_2_to_3))
  }
}


# Apply spline transformation (only for protein variables)
prepare_lasso_data <- function(data, predictors, spline_df = 3) {
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


# SHAP analysis function to replace variable importance
analyze_shap_values <- function(models, test_data, top_n = 50) {
  
  library(fastshap)
  
  shap_results <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Prepare prediction function for SHAP
    if (model_name == "xgb") {
      # For XGBoost model
      pfun <- function(object, newdata) {
        predict(object, as.matrix(newdata))
      }
    } else {
      # For other models (svm, rf, etc.)
      pfun <- function(object, newdata) {
        predict(object, newdata)
      }
    }
    
    # Calculate SHAP values
    shap_values <- fastshap::explain(
      model, 
      X = test_data[, -which(names(test_data) == "target")],
      pred_wrapper = pfun,
      nsim = 100  # Number of Monte Carlo simulations
    )
    
    # Calculate both absolute and signed SHAP values
    mean_abs_shap <- colMeans(abs(shap_values))  # For magnitude/importance
    mean_signed_shap <- colMeans(shap_values)    # For direction/impact
    
    # Create importance data frame with both magnitude and direction
    shap_importance <- data.frame(
      Feature = names(mean_abs_shap),
      Mean_Abs_SHAP = as.numeric(mean_abs_shap),
      Mean_Signed_SHAP = as.numeric(mean_signed_shap),
      Impact_Direction = ifelse(mean_signed_shap > 0, "Positive", "Negative")
    )
    
    # Sort by absolute importance (magnitude)
    shap_importance <- shap_importance[order(shap_importance$Mean_Abs_SHAP, decreasing = TRUE), ]
    
    shap_results[[model_name]] <- list(
      shap_values = shap_values,
      importance = head(shap_importance, top_n),
      model_type = model_name,
      raw_shap_matrix = shap_values  # Store raw values for detailed analysis
    )
  }
  
  # Combine SHAP importance across all models with direction information
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
  
  # Sort by absolute importance
  combined_importance <- combined_importance[order(combined_importance$Mean_Abs_SHAP, decreasing = TRUE), ]
  
  return(list(
    by_model = shap_results,
    average = combined_importance,
    top_features = head(combined_importance, top_n)
  ))
}


# Feature selection function for binary classification
feature_selection <- function(X, y, top_k = 40, final_k = 15) {
  # Calculate effect sizes using Cohen's d
  effect_sizes <- apply(X, 2, function(x) {
    abs(mean(x[y == 1]) - mean(x[y == 0])) / sd(x)
  })
  
  # Select top features based on effect size
  top_features <- names(sort(effect_sizes, decreasing = TRUE)[1:top_k])
  
  # Recursive feature elimination with random forest
  set.seed(123)
  ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
  rfe_result <- rfe(X[, top_features], y, 
                    sizes = c(10, 20, 30), 
                    rfeControl = ctrl)
  
  # LASSO regularization for final selection
  cv_fit <- cv.glmnet(as.matrix(X[, rfe_result$optVariables]), y, 
                      family = "binomial", alpha = 1)
  coefs <- coef(cv_fit, s = "lambda.min")
  selected_vars <- rownames(coefs)[which(coefs != 0)][-1]  # Remove intercept
  
  # Return final selected features
  if (length(selected_vars) > final_k) {
    return(head(selected_vars, final_k))
  } else {
    return(selected_vars)
  }
}


# LASSO regression for feature selection (continuous outcomes)
perform_lasso_selection <- function(train_data, alpha = 1) {
  
  x <- as.matrix(train_data[, -1])
  y <- train_data$target
  
  # Cross-validated LASSO
  cv_lasso <- cv.glmnet(x, y, 
                        family = "gaussian",  
                        alpha = alpha,
                        nfolds = 10,
                        type.measure = "mse") 
  
  best_lambda <- cv_lasso$lambda.1se
  
  # Fit final LASSO model
  lasso_model <- glmnet(x, y, 
                        family = "gaussian", 
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

# Compare regression model performance
compare_regression_models <- function(models, test_data) {
  performance_df <- data.frame()
  predictions_list <- list()
  
  for (model_name in names(models)) {
    # Generate predictions
    raw_predictions <- predict(models[[model_name]], test_data)
    predictions <- pmax(raw_predictions, 0.01)  # Set minimum value to 0.01
    true_values <- test_data$target
    
    # Calculate performance metrics
    rmse <- rmse(true_values, predictions)
    mae <- mae(true_values, predictions)
    r_squared <- cor(true_values, predictions)^2
    
    negative_ratio <- sum(predictions <= 0) / length(predictions)
    
    performance_df <- rbind(performance_df, data.frame(
      Model = model_name,
      RMSE = rmse,
      MAE = mae,
      R_squared = r_squared,
      Negative_Ratio = negative_ratio
    ))
    predictions_list[[model_name]] <- predictions
  }
  
  performance_df <- performance_df[order(-performance_df$R_squared), ]
  return(list(
    performance = performance_df,
    predictions = predictions_list
  ))
}

# Train multiple regression models
train_regression_models <- function(train_data) {
  
  models <- list()
  
  # Ensure positive target values
  if(any(train_data$target <= 0)) {
    train_data$target <- train_data$target - min(train_data$target) + 0.001
  }
  
  # Configure cross-validation
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    savePredictions = "final"
  )
  
  # 6.1 Generalized Linear Model
  cat("Training Generalized Linear Model...\n")
  
  models$glm <- train(
    target ~ .,
    data = train_data,
    trControl = ctrl,
    method = "glm",
    family = gaussian(link = "identity"),
    metric = "RMSE"
  )
  
  # 6.2 Random Forest
  cat("Training Random Forest Model...\n")
  models$rf <- train(
    target ~ .,
    data = train_data,
    method = "rf",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 3
  )
  
  # 6.3 XGBoost
  cat("Training XGBoost Model...\n")
  models$xgb <- train(
    target ~ .,
    data = train_data,
    method = "xgbTree",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 3,
    tuneGrid = expand.grid(
      nrounds = c(100, 150),
      max_depth = c(3, 4),
      eta = c(0.1, 0.3),
      gamma = 0,
      colsample_bytree = 0.8,
      min_child_weight = 1,
      subsample = 0.8
    )
  )
  
  # 6.4 Elastic Net
  cat("Training Elastic Net Model...\n")
  models$glmnet <- train(
    target ~ .,
    data = train_data,
    method = "glmnet",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5
  )
  
  # 6.5 Support Vector Machine
  cat("Training Support Vector Machine Model...\n")
  models$svm <- train(
    target ~ .,
    data = train_data,
    method = "svmRadial",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5
  )
  
  cat("6. Training Neural Network...\n")
  models$nnet <- train(
    target ~ .,
    data = train_data,
    method = "nnet",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5,
    trace = FALSE,
    linout = TRUE,
        tuneGrid = expand.grid(
      size = c(1,3,5,7,9),
      decay = c(0.0001,0.001,0.01,0.1,0.5)
    ),
    maxit = 200
  )
  
  # 7. Additional: Gradient Boosting Machine
  cat("Training Gradient Boosting Machine...\n")
  models$gbm <- train(
    target ~ .,
    data = train_data,
    method = "gbm",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 5,
    verbose = FALSE)
  
  return(models)
}

# Create ensemble model
create_ensemble <- function(predictions_list, true_values, method = "weighted") {
  
  if (method == "weighted") {
    # Weighted average based on R-squared performance
    weights <- comparison_results$performance$R_squared
    weights <- weights / sum(weights)  
    
    ensemble_pred <- matrix(0, nrow = length(predictions_list[[1]]))
    for (i in 1:length(predictions_list)) {
      ensemble_pred <- ensemble_pred + weights[i] * predictions_list[[i]]
    }
    
  } else if (method == "simple_average") {
    # Simple average of all predictions
    all_preds <- do.call(cbind, predictions_list)
    ensemble_pred <- rowMeans(all_preds)
  }
  
  # Ensure positive predictions
  ensemble_pred <- pmax(ensemble_pred, 0.01)
  
  # Calculate ensemble performance
  ensemble_rmse <- rmse(true_values, ensemble_pred)
  ensemble_mae <- mae(true_values, ensemble_pred)
  ensemble_r2 <- cor(true_values, ensemble_pred)^2
  
  results <- list(
    predictions = ensemble_pred,
    performance = data.frame(
      Model = "Ensemble",
      RMSE = ensemble_rmse,
      MAE = ensemble_mae,
      R_squared = ensemble_r2,
      Negative_Ratio = 0
    )
  )
  return(results)
}


# Build protein risk score using coefficients
build_protein_score <- function(data, coefficients) {
  
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

# Evaluate protein score performance
evaluate_protein_score <- function(protein_score, true_values, dataset_name = "") {
  
  # Calculate evaluation metrics
  rmse_val <- rmse(true_values, protein_score)
  mae_val <- mae(true_values, protein_score)
  r_squared <- cor(true_values, protein_score)^2
  correlation <- cor(true_values, protein_score)
  
  # Prepare data for visualization
  plot_data <- data.frame(
    True_Value = true_values,
    Predicted_Score = protein_score
  )
  
  results <- list(
    rmse = rmse_val,
    mae = mae_val,
    r_squared = r_squared,
    correlation = correlation,
    plot_data = plot_data
  )
  return(results)
}

# Create regression diagnostic plots
create_regression_plots <- function(train_results, test_results) {
  
  # Prepare training data for plotting
  train_plot_data <- train_results$plot_data
  train_plot_data$Dataset <- paste0("Training data (R² = ", round(train_results$r_squared, 3), ")")
  
  # Prepare test data for plotting
  test_plot_data <- test_results$plot_data
  test_plot_data$Dataset <- paste0("Testing data (R² = ", round(test_results$r_squared, 3), ")")
  
  # Combine datasets
  all_plot_data <- rbind(train_plot_data, test_plot_data)
  
  # Create scatter plot
  scatter_plot <- ggplot(all_plot_data, aes(x = True_Value, y = Predicted_Score, color = Dataset)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Training Set (R² = 0.5)" = "blue", 
                                  "Test Set (R² = 0.5)" = "red")) +
    labs(title = "Protein Score Prediction Performance",
         subtitle = "Actual Values vs Predicted Values",
         x = "Actual Fat Deposition Level",
         y = "Predicted Protein Score") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Create residual plot
  all_plot_data$Residuals <- all_plot_data$True_Value - all_plot_data$Predicted_Score
  residual_plot <- ggplot(all_plot_data, aes(x = Predicted_Score, y = Residuals, color = Dataset)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Training Set  (R² = 0.5)" = "blue", 
                                  "Test Set (R² = 0.5)" = "red")) +
    labs(title = "Residual Analysis",
         x = "Predicted Values",
         y = "Residuals") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(
    scatter_plot = scatter_plot,
    residual_plot = residual_plot
  ))
}

# Evaluate classification model performance
evaluate_model <- function(model, X_test, y_test) {
  if (model$type == "glmnet") {
    preds <- predict(model$model, newx = as.matrix(X_test[, model$features]), 
                     type = "response", s = "lambda.min")
  } else {
    preds <- predict(model, as.matrix(X_test[, model$features]))
  }
  
  # Calculate ROC and AUC
  roc_obj <- roc(y_test, as.numeric(preds))
  auc_value <- auc(roc_obj)
  
  # Calculate confusion matrix
  pred_class <- ifelse(preds > 0.5, 1, 0)
  conf_matrix <- confusionMatrix(factor(pred_class), factor(y_test))
  
  return(list(
    auc = auc_value,
    confusion = conf_matrix,
    roc = roc_obj
  ))
}


# Alternative approach: If ensemble models are still to be considered, model interpretability methods can be used
analyze_ensemble_importance <- function(models_list, model_index, test_data, 
                                        top_n = 20) {
  # Get ensemble model weights
  ensemble_weights <- models_list[[model_index]]$ensemble_weights
  
  if (is.null(ensemble_weights)) {
    # If weights are not stored, calculate performance-based weights
    model_perf <- models_list[[model_index]]$model_performance
    base_models <- model_perf[!grepl("Ensemble", model_perf$Model), ]
    
    # Calculate weights based on R²
    weights <- base_models$R_squared / sum(base_models$R_squared)
    names(weights) <- base_models$Model
  } else {
    weights <- ensemble_weights
  }
  
  # Calculate weighted feature importance
  weighted_importance <- list()
  
  for (model_name in names(weights)) {
    model_idx <- which(models_list[[model_index]]$model_performance$Model == model_name)
    if (length(model_idx) > 0) {
      # Get SHAP values or feature importance for this model
      model_shap <- models_list[[model_index]]$shap_analysis[[model_name]]
      
      if (!is.null(model_shap)) {
        # Calculate weighted importance
        for (feature in names(model_shap$feature_importance)) {
          if (is.null(weighted_importance[[feature]])) {
            weighted_importance[[feature]] <- 0
          }
          weighted_importance[[feature]] <- weighted_importance[[feature]] + 
            model_shap$feature_importance[[feature]] * weights[model_name]
        }
      }
    }
  }
}


# Complete main loop SHAP processing logic
process_shap_analysis <- function(all_models_results, tar_var) {
  shap_results <- list()
  
  for (i in 1:length(all_models_results)) {
    model_type <- all_models_results[[i]]$model_type
    model_desc <- all_models_results[[i]]$model_description
    
    # Check if the best model is an ensemble model
    best_perf <- all_models_results[[i]]$model_performance[which.max(all_models_results[[i]]$model_performance$R_squared), ]
    is_ensemble <- grepl("Ensemble", best_perf$Model, ignore.case = TRUE)
    
    if (is_ensemble) {
      cat("Model", model_type, "(", model_desc, ") best model is an ensemble model\n")
      
      # Use ensemble model interpretability analysis
      ensemble_analysis <- analyze_ensemble_importance(
        models_list = all_models_results,
        model_index = i,
        test_data = test_data,  # Need to pass test data
        top_n = 20
      )
      
      shap_results[[i]] <- list(
        type = "ensemble_analysis",
        plot = ensemble_analysis$plot,
        data = ensemble_analysis
      )
    } else {
      # Use regular SHAP analysis
      shap_plot <- plot_shap_summary(
        shap_results = all_models_results[[i]]$shap_analysis,
        model_name = paste0(model_desc, " (", tar_var, ")"),
        top_n = 20
      )
      
      shap_results[[i]] <- list(
        type = "shap_analysis",
        plot = shap_plot,
        best_model = best_perf$Model
      )
    }
  }
  return(shap_results)
}



# 2. Re-run LASSO on fixed feature set (without feature selection)
lasso_fixed_features <- function(data, features, alpha = 1, nfolds = 10) {
  
  # Ensure only pre-selected features are used
  X_subset <- data[, features]
  
  # Use cross-validation to select lambda, but set a lower minimum lambda
  cv_lasso <- cv.glmnet(
    x = as.matrix(X_subset),
    y = data$target,
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    type.measure = "mse",
    lambda.min.ratio = 0.001  # Set lower minimum lambda to reduce shrinkage
  )
  
  # Use lambda.min (instead of lambda.1se) to retain more coefficients
  coef_matrix <- as.matrix(coef(cv_lasso, s = "lambda.min"))
  coefficients <- coef_matrix[-1, 1]  # Remove intercept
  names(coefficients) <- features
  
  return(coefficients)
}


# Compare performance of three models
compare_three_models <- function(models_list) {
  performance_df <- data.frame()
  
  for (i in 1:length(models_list)) {
    model_perf <- models_list[[i]]$model_performance
    best_perf <- model_perf[which.max(model_perf$R_squared), ]
    
    performance_df <- rbind(performance_df, data.frame(
      Model_Type = models_list[[i]]$model_type,
      Model_Description = models_list[[i]]$model_description,
      Best_Model = best_perf$Model,
      R_squared = best_perf$R_squared,
      RMSE = best_perf$RMSE,
      MAE = best_perf$MAE,
      Feature_Count = models_list[[i]]$feature_count,
      Proteins_Count = length(models_list[[i]]$selected_proteins)
    ))
  }
  
  return(performance_df)
}


# Plot model comparison
plot_model_comparison <- function(performance_df) {
  p <- ggplot(performance_df, aes(x = Model_Description, y = R_squared, fill = Model_Description)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sprintf("R² = %.3f\nFeatures = %d", R_squared, Feature_Count)), 
              vjust = -0.3, size = 3) +
    scale_fill_jco() +
    labs(title = "Hierarchical Model Performance Comparison",
         x = "Model Type", y = "R-squared") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  return(p)
}

# Find the best model
find_best_model <- function(models_list) {
  best_r2 <- -Inf
  best_model_index <- 1
  best_model_name <- ""
  
  for (i in 1:length(models_list)) {
    model_perf <- models_list[[i]]$model_performance
    current_best_r2 <- max(model_perf$R_squared, na.rm = TRUE)
    
    if (current_best_r2 > best_r2) {
      best_r2 <- current_best_r2
      best_model_index <- i
      best_model_name <- model_perf$Model[which.max(model_perf$R_squared)]
    }
  }
  
  return(list(
    best_model_index = best_model_index,
    best_model_name = best_model_name,
    best_r_squared = best_r2
  ))
}

# Variance decomposition analysis
analyze_variance_decomposition <- function(models_list) {
  # Calculate additional variance explained by each hierarchical level
  base_r2 <- max(models_list[[1]]$model_performance$R_squared, na.rm = TRUE)
  basic_r2 <- max(models_list[[2]]$model_performance$R_squared, na.rm = TRUE)
  full_r2 <- max(models_list[[3]]$model_performance$R_squared, na.rm = TRUE)
  
  variance_breakdown <- data.frame(
    Component = c("Proteins", "Basic Information", "Body Composition", "Total Variance"),
    R_squared = c(base_r2, basic_r2 - base_r2, full_r2 - basic_r2, full_r2),
    Proportion = c(base_r2/full_r2, (basic_r2 - base_r2)/full_r2, 
                   (full_r2 - basic_r2)/full_r2, 1)
  )
  
  return(variance_breakdown)
}

# Print hierarchical model summary
print_hierarchical_model_summary <- function(models_list, target_var) {
  cat("\n")
  cat("===========================================\n")
  cat("Target Variable:", target_var, "\n")
  cat("===========================================\n")
  
  # Protein selection summary
  cat("Protein Selection Results:\n")
  cat("- Number of selected proteins:", length(models_list[[1]]$selected_proteins), "\n")
  cat("- Selected proteins:", paste(models_list[[1]]$selected_proteins, collapse = ", "), "\n")
  
  # Model performance summary
  cat("\nModel Performance Comparison:\n")
  for (i in 1:length(models_list)) {
    best_perf <- models_list[[i]]$model_performance[which.max(models_list[[i]]$model_performance$R_squared), ]
    cat(sprintf("- %s: R² = %.3f, RMSE = %.3f, Feature Count = %d\n",
                models_list[[i]]$model_description,
                best_perf$R_squared,
                best_perf$RMSE,
                models_list[[i]]$feature_count))
  }
  
  # Variance decomposition
  var_decomp <- analyze_variance_decomposition(models_list)
  cat("\nVariance Decomposition:\n")
  for (j in 1:3) {
    cat(sprintf("- %s: Explained Variance = %.3f (%.1f%%)\n",
                var_decomp$Component[j],
                var_decomp$R_squared[j],
                var_decomp$Proportion[j] * 100))
  }
  cat("===========================================\n\n")
}
