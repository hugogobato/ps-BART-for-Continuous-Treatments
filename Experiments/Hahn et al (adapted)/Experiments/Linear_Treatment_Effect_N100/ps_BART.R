number_experiments <- 100
number_of_points <- 250
data<-read.csv("dataset_1.csv")
num_samples <- length(data$y)

# Function to calculate the empirical CDF and map predictions to percentiles
map_to_percentiles <- function(observed, predictions) {
  # Calculate empirical CDF for observed data
  ecdf_function <- ecdf(observed)
  
  # Map predictions to their corresponding percentiles
  percentiles <- ecdf_function(predictions)
  
  return(percentiles)
}

# Function to calculate the loc part of Y
calculate_loc_Y <- function(X1, X2, X3, X4, X5, D) {
  
  mu_func <- function(X4, X1, X3) { -6 + X4 + 6 * abs(X3 - 1) }
  tau_func <- function(X2, X5) { 1 + 2 * X2 * X5 }
  
  mu <- mu_func(X4, X1, X3)
  tau <- tau_func(X2, X5)
  
  outcome <- mu + D * tau

  return(outcome)
}


library(dbarts)
n_posterior_simulations <- 1000

empirical_ATE_RMSE <-numeric(number_experiments)
empirical_ATE_MAE <-numeric(number_experiments)
empirical_ATE_MAPE <-numeric(number_experiments)
empirical_ATE_len <-numeric(number_experiments)
empirical_ATE_cover <-numeric(number_experiments)

theoretical_CATE_RMSE <- matrix(nrow = num_samples, ncol = number_experiments)
theoretical_CATE_MAE <-matrix(nrow = num_samples, ncol = number_experiments)
theoretical_CATE_MAPE <-matrix(nrow = num_samples, ncol = number_experiments)
theoretical_CATE_len <-matrix(nrow = num_samples, ncol = number_experiments)
theoretical_CATE_cover <-matrix(nrow = num_samples, ncol = number_experiments)

for (i in 1:number_experiments) {

file_name <- paste0("dataset_", i, ".csv")
data<-read.csv(file_name)
Y <- data$y
X <- data[, !(names(data) %in% c("y", "tau", "mu", "D"))]
D <- data$D

max_D <- max(D)

dbfit <- bart(X,D,X,verbose=F)
mean_predictions <- colMeans(dbfit$yhat.test)

# Map the mean predictions to percentiles using the empirical CDF of D
pi_hats <- map_to_percentiles(D, mean_predictions)

X <- cbind(X,pi_hats)

# Calculate the column means of X
column_means <- colMeans(X)

# Create a matrix where all values in each column are equal to the column mean
X_mean <- matrix(rep(column_means, each = nrow(X)), nrow = nrow(X), byrow = FALSE)

bart_model <- bart(cbind(X,D),Y,cbind(X_mean,D),verbose=F, keeptrees = TRUE)

predicted_y <- colMeans(bart_model$yhat.test)

real_ATE <- calculate_loc_Y(X_mean[,1], X_mean[,2],X_mean[,3], X_mean[,4],X_mean[,5], D)

empirical_ATE_RMSE[i] <- sqrt(mean((real_ATE - predicted_y)^2))
empirical_ATE_MAE[i] <- mean(abs(real_ATE - predicted_y))
empirical_ATE_MAPE[i] <- mean(abs(real_ATE - predicted_y)/abs(real_ATE))

# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
empirical_ATE_lower_bounds <- apply(bart_model$yhat.test, 2, function(x) quantile(x, probs = 0.025))
empirical_ATE_upper_bounds <- apply(bart_model$yhat.test, 2, function(x) quantile(x, probs = 0.975))

empirical_ATE_len[i] <- mean(abs(empirical_ATE_upper_bounds-empirical_ATE_lower_bounds))

cover_vector <- (real_ATE >= empirical_ATE_lower_bounds) & (real_ATE <= empirical_ATE_upper_bounds)

# Calculate the number of times tau is within the CI
number_covered <- sum(cover_vector)

empirical_ATE_cover[i] <- number_covered / length(real_ATE)

for (j in 1:num_samples){
column_means[1] <- X$x1[j]
column_means[2] <- X$x2[j]
column_means[3] <- X$x3[j]
column_means[4] <- X$x4[j]
column_means[5] <- X$x5[j]
column_means[6] <- X$pi_hats[j]
D_range <- seq(0, max_D, length.out = number_of_points)

X_mean <- matrix(rep(column_means, each = number_of_points), nrow = number_of_points, byrow = FALSE)
loc_Y_D <- calculate_loc_Y(X_mean[,1], X_mean[,2],X_mean[,3], X_mean[,4],X_mean[,5], D_range)

new_data <- cbind(X_mean, D_range)
all_estimated_theoretical_y <- predict(bart_model, newdata = new_data,)
estimated_theoretical_y <- colMeans(all_estimated_theoretical_y)

theoretical_CATE_RMSE[j,i] <- sqrt(mean((loc_Y_D - estimated_theoretical_y)^2))
theoretical_CATE_MAE[j,i] <- mean(abs(loc_Y_D - estimated_theoretical_y))
theoretical_CATE_MAPE[j,i] <- mean(abs(loc_Y_D - estimated_theoretical_y)/abs(loc_Y_D))

# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
theoretical_CATE_lower_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.025))
theoretical_CATE_upper_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.975))

theoretical_CATE_len[j,i] <- mean(abs(theoretical_CATE_upper_bounds-theoretical_CATE_lower_bounds))

cover_vector <- (loc_Y_D >= theoretical_CATE_lower_bounds) & (loc_Y_D <= theoretical_CATE_upper_bounds)

# Calculate the number of times tau is within the CI
number_covered <- sum(cover_vector)

theoretical_CATE_cover[j,i] <- number_covered / length(loc_Y_D)

}

}

cat("The empirical ATE RMSE is:", mean(empirical_ATE_RMSE), "\n")
cat("The empirical ATE MAE is:", mean(empirical_ATE_MAE), "\n")
cat("The empirical ATE MAPE is:", mean(empirical_ATE_MAPE), "\n")
cat("The empirical ATE Length is:", mean(empirical_ATE_len), "\n")
cat("The empirical ATE Cover is:", mean(empirical_ATE_cover), "\n")

cat("The theoretical CATE RMSE is:", mean(colMeans(theoretical_CATE_RMSE)), "\n")
cat("The theoretical CATE MAE is:", mean(colMeans(theoretical_CATE_MAE)), "\n")
cat("The theoretical CATE MAPE is:", mean(colMeans(theoretical_CATE_MAPE)), "\n")
cat("The theoretical CATE Length is:", mean(colMeans(theoretical_CATE_len)), "\n")
cat("The theoretical CATE Cover is:", mean(colMeans(theoretical_CATE_cover)), "\n")

df <- data.frame(
  empirical_ATE_RMSE = empirical_ATE_RMSE,
  empirical_ATE_MAE = empirical_ATE_MAE,
  empirical_ATE_MAPE = empirical_ATE_MAPE,
  empirical_ATE_len = empirical_ATE_len,
  empirical_ATE_cover = empirical_ATE_cover,
  theoretical_CATE_RMSE=colMeans(theoretical_CATE_RMSE),
  theoretical_CATE_MAE=colMeans(theoretical_CATE_MAE),
  theoretical_CATE_MAPE=colMeans(theoretical_CATE_MAPE),
  theoretical_CATE_len=colMeans(theoretical_CATE_len),
  theoretical_CATE_cover=colMeans(theoretical_CATE_cover)
)

library(openxlsx)
file_path <- "ps_BART.xlsx"
write.xlsx(df, file_path, row.names = FALSE)


# Create a new workbook
wb <- createWorkbook()

# Add sheets with matrices
addWorksheet(wb, "CATE_RMSE")
writeData(wb, "CATE_RMSE", theoretical_CATE_RMSE)

addWorksheet(wb, "CATE_MAE")
writeData(wb, "CATE_MAE", theoretical_CATE_MAE)

addWorksheet(wb, "CATE_MAPE")
writeData(wb, "CATE_MAPE", theoretical_CATE_MAPE)

addWorksheet(wb, "CATE_len")
writeData(wb, "CATE_len", theoretical_CATE_len)

addWorksheet(wb, "CATE_cover")
writeData(wb, "CATE_cover", theoretical_CATE_cover)

# Save the workbook as an Excel file
saveWorkbook(wb, "ps_BART_CATE_full.xlsx", overwrite = TRUE)