library(dplyr)

load(file = "Process_Data/Fat_Dis/Model_Data.RData")

summarize_transformation_guidelines <- function(target_vector) {
  
  # Calculate comprehensive statistics
  stats <- list(
    n = length(target_vector),
    mean = mean(target_vector, na.rm = TRUE),
    median = median(target_vector, na.rm = TRUE),
    sd = sd(target_vector, na.rm = TRUE),
    skewness = moments::skewness(target_vector, na.rm = TRUE),
    cv = sd(target_vector, na.rm = TRUE) / mean(target_vector, na.rm = TRUE),
    zero_count = sum(target_vector == 0, na.rm = TRUE),
    negative_count = sum(target_vector < 0, na.rm = TRUE),
    range = range(target_vector, na.rm = TRUE),
    iqr = IQR(target_vector, na.rm = TRUE),
    q1 = quantile(target_vector, 0.25, na.rm = TRUE),
    q3 = quantile(target_vector, 0.75, na.rm = TRUE)
  )
  
  # Calculate additional derived metrics
  stats$range_ratio <- stats$range[2] / stats$range[1]
  stats$mean_median_ratio <- stats$mean / stats$median
  
  cat("=== TRANSFORMATION GUIDELINES SUMMARY ===\n\n")
  
  # Basic distribution statistics
  cat("BASIC STATISTICS:\n")
  cat("• n =", stats$n, "observations\n")
  cat("• Mean =", round(stats$mean, 3), "| Median =", round(stats$median, 3), "\n")
  cat("• SD =", round(stats$sd, 3), "| IQR =", round(stats$iqr, 3), "\n")
  cat("• Range: [", round(stats$range[1], 3), ",", round(stats$range[2], 3), "]\n")
  cat("• Quartiles: Q1 =", round(stats$q1, 3), "| Q3 =", round(stats$q3, 3), "\n\n")
  
  # Distribution shape assessment
  cat("DISTRIBUTION SHAPE:\n")
  
  # Skewness guidance with more granular categories
  if (abs(stats$skewness) > 2) {
    direction <- ifelse(stats$skewness > 0, "right", "left")
    cat("• Skewness =", round(stats$skewness, 3), ": STRONG", direction, "skew → STRONGLY RECOMMEND transformation\n")
  } else if (abs(stats$skewness) > 1) {
    direction <- ifelse(stats$skewness > 0, "right", "left")
    cat("• Skewness =", round(stats$skewness, 3), ": Moderate", direction, "skew → RECOMMEND transformation\n") 
  } else if (abs(stats$skewness) > 0.5) {
    direction <- ifelse(stats$skewness > 0, "right", "left")
    cat("• Skewness =", round(stats$skewness, 3), ": Mild", direction, "skew → CONSIDER transformation\n")
  } else {
    cat("• Skewness =", round(stats$skewness, 3), ": Approximately symmetric → Usually no transformation needed\n")
  }
  
  # Mean-median comparison for additional skewness insight
  if (stats$mean_median_ratio > 1.1) {
    cat("• Mean > Median suggests right skew (ratio =", round(stats$mean_median_ratio, 3), ")\n")
  } else if (stats$mean_median_ratio < 0.9) {
    cat("• Mean < Median suggests left skew (ratio =", round(stats$mean_median_ratio, 3), ")\n")
  } else {
    cat("• Mean ≈ Median suggests symmetric distribution (ratio =", round(stats$mean_median_ratio, 3), ")\n")
  }
  
  cat("\nVARIANCE CHARACTERISTICS:\n")
  
  # Coefficient of variation guidance
  if (stats$cv > 1) {
    cat("• High CV =", round(stats$cv, 3), ": Variance increases with mean → Transformation may improve stability\n")
  } else if (stats$cv > 0.5) {
    cat("• Moderate CV =", round(stats$cv, 3), ": Some variance-mean relationship present\n")
  } else {
    cat("• Low CV =", round(stats$cv, 3), ": Relatively constant variance\n")
  }
  
  cat("\nDATA COMPOSITION:\n")
  
  # Zero and negative values guidance
  if (stats$zero_count > 0) {
    zero_percent <- round(100 * stats$zero_count / stats$n, 1)
    cat("• Contains", stats$zero_count, "zero values (", zero_percent, "%): Consider log1p or sqrt transformation\n")
  } else {
    cat("• No zero values: Standard log transformation is feasible\n")
  }
  
  if (stats$negative_count > 0) {
    neg_percent <- round(100 * stats$negative_count / stats$n, 1)
    cat("• Contains", stats$negative_count, "negative values (", neg_percent, "%): Consider Box-Cox (with offset) or keep original\n")
  } else {
    cat("• No negative values\n")
  }
  
  cat("\nDYNAMIC RANGE ASSESSMENT:\n")
  
  # Range guidance
  if (stats$range_ratio > 1000 && stats$range[1] > 0) {
    cat("• Large dynamic range (ratio =", round(stats$range_ratio, 1), "): Log transformation is USUALLY BENEFICIAL\n")
  } else if (stats$range_ratio > 100) {
    cat("• Moderate dynamic range (ratio =", round(stats$range_ratio, 1), "): Log transformation may help\n")
  } else {
    cat("• Limited dynamic range (ratio =", round(stats$range_ratio, 1), "): Transformation less critical\n")
  }
  
  # IQR to range ratio for outlier assessment
  iqr_range_ratio <- stats$iqr / (stats$range[2] - stats$range[1])
  if (iqr_range_ratio < 0.3) {
    cat("• Low IQR/Range ratio (", round(iqr_range_ratio, 3), ") suggests potential outliers\n")
  }
  
  cat("\nRECOMMENDED TRANSFORMATIONS (in order of preference):\n")
  
  # Generate transformation recommendations based on characteristics
  if (stats$zero_count > 0 && stats$negative_count == 0) {
    cat("1. log1p(x) = log(x + 1)  [handles zeros]\n")
    cat("2. sqrt(x)  [handles zeros, milder effect]\n")
  } else if (stats$negative_count > 0) {
    cat("1. Box-Cox with appropriate offset\n") 
    cat("2. Keep original scale or consider alternative models\n")
  } else if (stats$skewness > 1 && stats$range[1] > 0) {
    cat("1. log(x)  [for right-skewed positive data]\n")
    cat("2. sqrt(x)  [milder alternative]\n")
    cat("3. 1/x  [for severe right skew]\n")
  } else if (stats$skewness < -1) {
    cat("1. x² or x³  [for left-skewed data]\n")
    cat("2. Keep original scale\n")
  } else {
    cat("• No strong transformation needed based on current metrics\n")
  }
  
  return(invisible(stats))
}

# Apply to all four variables
cat("Analyzing Pancreas_PDFF_fat_fraction:\n")
guidelines_pancreas <- summarize_transformation_guidelines(dat.ana.u.1$Pancreas_PDFF_fat_fraction)
cat("\n")

cat("Analyzing Liver_PDFF_fat_fraction:\n")
guidelines_liver <- summarize_transformation_guidelines(dat.ana.u.1$Liver_PDFF_fat_fraction)
cat("\n")

cat("Analyzing Muscle_fat_infiltration:\n")
guidelines_muscle <- summarize_transformation_guidelines(dat.ana.u.1$Muscle_fat_infiltration)
cat("\n")

cat("Analyzing area_of_pericardial_fat:\n")
guidelines_pericardial <- summarize_transformation_guidelines(dat.ana.u.1$area_of_pericardial_fat)
cat("\n")

# Create a comprehensive summary table for easy comparison
variables_summary <- data.frame(
  Variable = c("Pancreas_PDFF_fat_fraction", "Liver_PDFF_fat_fraction", 
               "Muscle_fat_infiltration", "area_of_pericardial_fat"),
  n = c(guidelines_pancreas$n, guidelines_liver$n,
        guidelines_muscle$n, guidelines_pericardial$n),
  Mean = c(guidelines_pancreas$mean, guidelines_liver$mean,
           guidelines_muscle$mean, guidelines_pericardial$mean),
  Median = c(guidelines_pancreas$median, guidelines_liver$median,
             guidelines_muscle$median, guidelines_pericardial$median),
  SD = c(guidelines_pancreas$sd, guidelines_liver$sd,
         guidelines_muscle$sd, guidelines_pericardial$sd),
  Skewness = c(guidelines_pancreas$skewness, guidelines_liver$skewness,
               guidelines_muscle$skewness, guidelines_pericardial$skewness),
  CV = c(guidelines_pancreas$cv, guidelines_liver$cv,
         guidelines_muscle$cv, guidelines_pericardial$cv),
  Q1 = c(guidelines_pancreas$q1, guidelines_liver$q1,
         guidelines_muscle$q1, guidelines_pericardial$q1),
  Q3 = c(guidelines_pancreas$q3, guidelines_liver$q3,
         guidelines_muscle$q3, guidelines_pericardial$q3),
  IQR = c(guidelines_pancreas$iqr, guidelines_liver$iqr,
          guidelines_muscle$iqr, guidelines_pericardial$iqr),
  Min = c(guidelines_pancreas$range[1], guidelines_liver$range[1],
          guidelines_muscle$range[1], guidelines_pericardial$range[1]),
  Max = c(guidelines_pancreas$range[2], guidelines_liver$range[2],
          guidelines_muscle$range[2], guidelines_pericardial$range[2]),
  Range_Ratio = c(guidelines_pancreas$range_ratio, guidelines_liver$range_ratio,
                  guidelines_muscle$range_ratio, guidelines_pericardial$range_ratio),
  Mean_Median_Ratio = c(guidelines_pancreas$mean_median_ratio, guidelines_liver$mean_median_ratio,
                        guidelines_muscle$mean_median_ratio, guidelines_pericardial$mean_median_ratio),
  Zero_Count = c(guidelines_pancreas$zero_count, guidelines_liver$zero_count,
                 guidelines_muscle$zero_count, guidelines_pericardial$zero_count),
  Negative_Count = c(guidelines_pancreas$negative_count, guidelines_liver$negative_count,
                     guidelines_muscle$negative_count, guidelines_pericardial$negative_count)
)

# Format the summary table for better readability
formatted_summary <- variables_summary
numeric_columns <- sapply(formatted_summary, is.numeric)
formatted_summary[numeric_columns] <- round(formatted_summary[numeric_columns], 3)

# Display the formatted summary
print("=== COMPREHENSIVE VARIABLES SUMMARY ===")
print(formatted_summary)

# Add transformation recommendation based on skewness
formatted_summary$Transformation_Recommendation <- ifelse(
  abs(formatted_summary$Skewness) > 2, "Strongly Recommended",
  ifelse(abs(formatted_summary$Skewness) > 1, "Recommended",
         ifelse(abs(formatted_summary$Skewness) > 0.5, "Consider", "Not Needed"))
)

# Display final recommendations
cat("\n=== TRANSFORMATION RECOMMENDATIONS SUMMARY ===\n")
for(i in 1:nrow(formatted_summary)) {
  cat(sprintf("%s: %s (Skewness = %.3f)\n", 
              formatted_summary$Variable[i], 
              formatted_summary$Transformation_Recommendation[i],
              formatted_summary$Skewness[i]))
}
print("Summary Table for All Variables:")
print(formatted_summary)
write.csv(formatted_summary,file = "Process_Data/Fat_Dis/Fat_transformation.csv",row.names = F)

