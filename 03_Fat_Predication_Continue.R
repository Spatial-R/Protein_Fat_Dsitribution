
rm(list = ls())

library(dplyr)
library(stringi)
library(cowplot)
library(ggsci)
library(ggplot2)
library(tidyr)
library(readr)
library(VIM)
library(caret)
library(splines)
library(glmnet)
library(randomForest)
library(xgboost)
library(pROC)
library(ggpubr)
library(Metrics)
library(plotROC)

source("Codes/00_Basic_Functions.R")
source("Codes_Fat/00_ML_Continue_functions.R")
load("Process_Data/Fat_Dis/Fat_Pop.RData")
load("Process_Data/Fat_Dis/Pro_Impt.RData")

dat_body_fin <- filter(dat_body_fin_u, instance == 2)

tar_var <- c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction",
             "Muscle_fat_infiltration", "area_of_pericardial_fat")

tar_var_name <- c("age", "sex", "waist_circum", "BMI", tar_var)
dat.ana_2 <- dat.ana.1[, match(c("eid", tar_var_name), names(dat.ana.1))]

dat.ana_2 %>%
  inner_join(pro_int_2_2, by = c("eid")) %>%
  inner_join(dat_body_fin, by = c("eid")) -> dat.ana.u

ncol(pro_int_2_2) - 1

dat.ana.u <- mutate(dat.ana.u, sex = ifelse(sex == "Male", 1, 0))
dat.ana.u.1 <- fast_impute(dat.ana.u, method = "rf")

na.row <- data.frame(colSums(is.na(dat.ana.u.1)));names(na.row) <- "na_num"

dat.ana.u.1$instance <- NULL
protein_names <- names(dat.ana.u.1)[10:(ncol(dat.ana.u.1) - 27)]
basic_vars <- c("age", "sex", "waist_circum", "BMI")
body_comp_vars <- setdiff(names(dat_body_fin), c("eid", "instance", basic_vars))

###################  for the zero, we used 0.001   ##################
dat.ana.u.1 <- mutate(dat.ana.u.1, sex = as.character(sex))
dat.ana.u.1 <- filter(dat.ana.u.1,! Liver_PDFF_fat_fraction == 0)
dat.ana.u.1[,c(6:7)] <- apply(dat.ana.u.1[,c(6:7)],2,function(data){
  log(ifelse(data == 0,0.001,data))
})


# Modified main loop
for (i in c(1:length(tar_var))) {
  
  tar_col_tem <- match(tar_var[i], names(dat.ana.u.1))
  dat_tem <- dat.ana.u.1[, c(tar_col_tem, 2:5, 10:ncol(dat.ana.u.1))]
  names(dat_tem)[1] <- "target"
  
  set.seed(20251111)
  # Split training and test sets
  train_index <- createDataPartition(dat_tem$target, p = 0.7, list = FALSE)
  
  # Step 1: LASSO screening based on proteins only
  cat("=== Protein Screening Phase ===\n")
  cat("Target variable:", tar_var[i], "\n")
  
  # Prepare data containing only proteins
  protein_data <- dat_tem[, c("target", protein_names)]
  
  # Data preprocessing
  protein_ns_dat <- prepare_lasso_data(data = protein_data, 
                                       predictors = protein_names, spline_df = 3)
  protein_ns_dat1 <- cbind(protein_data["target"], protein_ns_dat$X)
  protein_processed <- preprocess_data(data = protein_ns_dat1)
  
  # LASSO feature selection
  lasso_results <- perform_lasso_selection(protein_processed)
  selected_proteins <- lasso_results$selected_features
  
  cat("Number of proteins selected by LASSO:", length(selected_proteins), "\n")
  
  # If no features are selected, use cross-validation to select top N important features
  if (length(selected_proteins) == 0) {
    cat("LASSO did not select any features, using cross-validation to select top 10 important features\n")
    full_coef <- as.matrix(coef(lasso_results$cv_model, s = lasso_results$cv_model$lambda.min))
    coef_df <- data.frame(
      Feature = rownames(full_coef),
      Coefficient = as.numeric(full_coef)
    )
    coef_df <- coef_df[coef_df$Feature != "(Intercept)", ]
    coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
    selected_proteins <- coef_df$Feature[1:min(10, nrow(coef_df))]
  }
  
  # Nonlinear effect analysis
  nonlinear_results <- analyze_nonlinear_effects(selected_proteins)
  target_proteins <- nonlinear_results$all_nonlinear
  cat("Proteins requiring spline transformation:", paste(target_proteins, collapse = ", "), "\n")
  
  # Store results for three models
  all_models_results <- list()
  
  for (model_type in 1:3) {
    
    cat("Processing target variable:", tar_var[i], "Model type:", model_type, "\n")
    
    # Prepare data according to model type
    if (model_type == 1) {
      # Model 1: Proteins only
      model_vars <- c("target", target_proteins)
      model_desc <- "Pure Protein Model"
    } else if (model_type == 2) {
      # Model 2: Proteins + Basic information
      model_vars <- c("target", target_proteins, basic_vars)
      model_desc <- "Protein + Basic Info Model"
    } else {
      # Model 3: Proteins + Basic information + Body composition
      model_vars <- c("target", target_proteins, basic_vars, body_comp_vars)
      model_desc <- "Full Model (Protein + Basic Info + Body Composition)"
    }
    
    # Select data subset
    dat_subset <- dat_tem[, match(model_vars, names(dat_tem))]
    
    # Data preprocessing
    data_processed <- preprocess_data(data = dat_subset)
    
    # Split training and test sets
    train_data <- data_processed[train_index, ]
    test_data <- data_processed[-train_index, ]
    
    # Train models
    models <- train_regression_models(train_data)
    
    # Model comparison
    comparison_results <- compare_regression_models(models, test_data)
    
    # Ensemble learning
    ensemble_results <- create_ensemble(
      comparison_results$predictions,
      test_data$target,
      method = "weighted"
    )
    
    # Variable importance analysis
    variable_importance <- analyze_variable_importance(models, top_n = 50)
    
    # SHAP analysis (replacing variable importance)
    shap_analysis <- analyze_shap_values(models, test_data, top_n = 50)
    
    # Protein score (only for Model 1)
    if (model_type == 1) {
      # Build protein score using fixed LASSO coefficients
      lasso_coef_fixed <- lasso_fixed_features(data = train_data, features = target_proteins)
      train_data$Protein_Score <- build_protein_score(train_data, lasso_coef_fixed)
      test_data$Protein_Score <- build_protein_score(test_data, lasso_coef_fixed)
      
      protein_score_train <- evaluate_protein_score(train_data$Protein_Score, train_data$target)
      protein_score_test <- evaluate_protein_score(test_data$Protein_Score, test_data$target)
    } else {
      protein_score_train <- NULL
      protein_score_test <- NULL
    }
    
    # Store current model results
    all_models_results[[model_type]] <- list(
      model_type = model_type,
      model_description = model_desc,
      selected_proteins = selected_proteins,
      models = models,
      origin_data = list(train_data, test_data),
      model_performance = rbind(comparison_results$performance, ensemble_results$performance),
      shap_analysis = shap_analysis,
      variable_importance = variable_importance,
      protein_score_performance = list(train = protein_score_train, test = protein_score_test),
      feature_count = ncol(dat_subset) - 1,
      variable_list = setdiff(model_vars, "target")
    )
  }
  
  # # Compare performance of three models
  # performance_comparison <- compare_three_models(all_models_results)
  # 
  # # Plot performance comparison
  # performance_plot <- plot_model_comparison(performance_comparison)

  # # Find best model and plot SHAP
  # best_model_info <- find_best_model(all_models_results)
  # 
  # shap_plot <- plot_shap_summary(
  #   shap_results = all_models_results[[best_model_info$best_model_index]]$shap_analysis,
  #   model_name = best_model_info$best_model_name, top_n = 20)
  # 
  # # Variance decomposition analysis
  # variance_decomposition <- analyze_variance_decomposition(all_models_results)

  # Save all results
  save(all_models_results,
       file = paste0("Process_Data/Fat_Dis/Continues/", tar_var[i], "_Hierarchical_Models.RData"))
  
  # Print summary report
  print_hierarchical_model_summary(all_models_results, tar_var[i])
}
