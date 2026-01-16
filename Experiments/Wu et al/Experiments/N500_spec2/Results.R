# Load necessary library
library(readxl)

# Read the Excel files
bcf_data <- read_excel("BCF.xlsx")
ps_bart_data <- read_excel("ps_BART.xlsx")

# Calculate means for BCF
bcf_means <- colMeans(bcf_data, na.rm = TRUE)

# Calculate means for ps-BART
ps_bart_means <- colMeans(ps_bart_data, na.rm = TRUE)

# Print results for BCF
cat("BCF:\n")
cat("The empirical ATE RMSE is:", bcf_means["empirical_ATE_RMSE"], "\n")
cat("The empirical ATE MAE is:", bcf_means["empirical_ATE_MAE"], "\n")
cat("The empirical ATE MAPE is:", bcf_means["empirical_ATE_MAPE"], "\n")
cat("The empirical ATE Length is:", bcf_means["empirical_ATE_len"], "\n")
cat("The empirical ATE Cover is:", bcf_means["empirical_ATE_cover"], "\n")
cat("The theoretical CATE RMSE is:", bcf_means["theoretical_CATE_RMSE"], "\n")
cat("The theoretical CATE MAE is:", bcf_means["theoretical_CATE_MAE"], "\n")
cat("The theoretical CATE MAPE is:", bcf_means["theoretical_CATE_MAPE"], "\n")
cat("The theoretical CATE Length is:", bcf_means["theoretical_CATE_len"], "\n")
cat("The theoretical CATE Cover is:", bcf_means["theoretical_CATE_cover"], "\n")

# Print results for ps-BART
cat("\nps-BART:\n")
cat("The empirical ATE RMSE is:", ps_bart_means["empirical_ATE_RMSE"], "\n")
cat("The empirical ATE MAE is:", ps_bart_means["empirical_ATE_MAE"], "\n")
cat("The empirical ATE MAPE is:", ps_bart_means["empirical_ATE_MAPE"], "\n")
cat("The empirical ATE Length is:", ps_bart_means["empirical_ATE_len"], "\n")
cat("The empirical ATE Cover is:", ps_bart_means["empirical_ATE_cover"], "\n")
cat("The theoretical CATE RMSE is:", ps_bart_means["theoretical_CATE_RMSE"], "\n")
cat("The theoretical CATE MAE is:", ps_bart_means["theoretical_CATE_MAE"], "\n")
cat("The theoretical CATE MAPE is:", ps_bart_means["theoretical_CATE_MAPE"], "\n")
cat("The theoretical CATE Length is:", ps_bart_means["theoretical_CATE_len"], "\n")
cat("The theoretical CATE Cover is:", ps_bart_means["theoretical_CATE_cover"], "\n")
