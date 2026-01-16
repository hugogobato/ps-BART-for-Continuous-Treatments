library(stochtree)
library(dbarts)

number_experiments <- 1
number_of_points <- 500
num_gfr <- 40
num_burnin <- 0
num_mcmc <- 1000

data<-read.csv("dataset_1.csv")
num_samples <- length(data$Y)

# Function to calculate the empirical CDF and map predictions to percentiles
map_to_percentiles <- function(observed, predictions) {
  # Calculate empirical CDF for observed data
  ecdf_function <- ecdf(observed)
  
  # Map predictions to their corresponding percentiles
  percentiles <- ecdf_function(predictions)
  
  return(percentiles)
}

# Function to calculate the loc part of Y
calculate_loc_Y <- function(X1, X2, D) {
  D + (X1 + X2) * exp(-D * (X1 + X2))
}

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
Y <- data$Y
X <- data[, !(names(data) %in% c("Y", "E_Y","D"))]
D <- data$D

max_D <- max(D)

dbfit <- bart(X,D,X,verbose=F)
mean_predictions <- colMeans(dbfit$yhat.test)

# Map the mean predictions to percentiles using the empirical CDF of D
pi_hats <- map_to_percentiles(D, mean_predictions)

X_full <- cbind(X,pi_hats)

# Calculate the column means of X
column_means <- colMeans(X_full)

# Create a matrix where all values in each column are equal to the column mean
X_mean <- matrix(rep(column_means, each = nrow(X)), nrow = nrow(X), byrow = FALSE)

dbfit <- bcf(
    X_train = X, Z_train = D, y_train = Y, , pi_train = pi_hats, num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc)

sample_inds <- (num_gfr + num_burnin + 1):num_mcmc

y_hat_test <- predict(dbfit, X_test= cbind(X_mean[,1], X_mean[,2]), Z_test = D, pi_test = X_mean[,3])

predicted_y <- rowMeans(y_hat_test$y_hat[,sample_inds])

real_ATE <- calculate_loc_Y(X_mean[,1], X_mean[,2], D)

empirical_ATE_RMSE[i] <- sqrt(mean((real_ATE - predicted_y)^2))
empirical_ATE_MAE[i] <- mean(abs(real_ATE - predicted_y))
empirical_ATE_MAPE[i] <- mean(abs(real_ATE - predicted_y)/abs(real_ATE))

# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
empirical_ATE_lower_bounds <- apply(y_hat_test$y_hat[,sample_inds], 1, function(x) quantile(x, probs = 0.025))
empirical_ATE_upper_bounds <- apply(y_hat_test$y_hat[,sample_inds], 1, function(x) quantile(x, probs = 0.975))

empirical_ATE_len[i] <- mean(abs(empirical_ATE_upper_bounds-empirical_ATE_lower_bounds))

cover_vector <- (real_ATE >= empirical_ATE_lower_bounds) & (real_ATE <= empirical_ATE_upper_bounds)

# Calculate the number of times tau is within the CI
number_covered <- sum(cover_vector)

empirical_ATE_cover[i] <- number_covered / length(real_ATE)

D_range <- seq(0, max_D, length.out = number_of_points)

loc_Y_D <- calculate_loc_Y(X_mean[,1], X_mean[,2], D_range)

# Define the y-range to cover all values
y_min <- min(c(predicted_y, empirical_ATE_lower_bounds, empirical_ATE_upper_bounds, loc_Y_D))
y_max <- max(c(predicted_y, empirical_ATE_lower_bounds, empirical_ATE_upper_bounds, loc_Y_D))

png("BCF_observed_D.png", width = 800, height = 600)
plot(D, predicted_y, 
     main = "Observed D (BCF)", 
     xlab = "D", 
     ylab = "Y",
     pch = 16,  # solid circles
     col = "blue",
     ylim = c(y_min, y_max))

# Add error bars
arrows(D, empirical_ATE_lower_bounds, D, empirical_ATE_upper_bounds, 
       angle = 90, code = 3, length = 0.05, col = "blue")

# Add the line plot of D_range and loc_Y_D
lines(D_range, loc_Y_D, type = "l", col = "red", lwd = 2)

# Add a legend to distinguish between the scatter plot and the line plot
legend("topright", legend = c("Estimated", "Actual"), 
       col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2))
dev.off()

X_mean <- matrix(rep(column_means, each = number_of_points), nrow = number_of_points, byrow = FALSE)
new_data <- cbind(X_mean, D_range)
predicted_y_CATE <- predict(dbfit, X_test= cbind(X_mean[,1], X_mean[,2]), Z_test = D_range,pi_test = X_mean[,3])
all_estimated_theoretical_y <- predicted_y_CATE$y_hat[,sample_inds]
estimated_theoretical_y <- rowMeans(all_estimated_theoretical_y)
# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
theoretical_CATE_lower_bounds <- apply(all_estimated_theoretical_y, 1, function(x) quantile(x, probs = 0.025))
theoretical_CATE_upper_bounds <- apply(all_estimated_theoretical_y, 1, function(x) quantile(x, probs = 0.975))

# Define the y-range to cover all values
y_min <- min(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
y_max <- max(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
png("BCF_ATE.png", width = 800, height = 600)
plot(D_range, estimated_theoretical_y, 
     main = "Estimated ATE (BCF)", 
     xlab = "D", 
     ylab = "Y",
     pch = 16,  # solid circles
     col = "blue",
     ylim = c(y_min, y_max))

# Add dotted lines for the lower and upper bounds
lines(D_range, theoretical_CATE_lower_bounds, type = "l", col = "blue", lwd = 1, lty = 2)
lines(D_range, theoretical_CATE_upper_bounds, type = "l", col = "blue", lwd = 1, lty = 2)

# Add the line plot of D_range and loc_Y_D
lines(D_range, loc_Y_D, type = "l", col = "red", lwd = 2)

# Add a legend to distinguish between the scatter plot and the line plot
legend("topright", legend = c("Estimated", "Actual", "Lower Bound", "Upper Bound"), 
       col = c("blue", "red", "blue", "blue"), pch = c(16, NA, NA, NA), 
       lty = c(NA, 1, 2, 2), lwd = c(NA, 2, 1, 1))
dev.off()

for (j in 1:1){
column_means[1] <- X_full$X1[j]
column_means[2] <- X_full$X2[j]
column_means[3] <- X_full$pi_hats[j]
D_range <- seq(0, max_D, length.out = number_of_points)

X_mean <- matrix(rep(column_means, each = number_of_points), nrow = number_of_points, byrow = FALSE)
loc_Y_D <- calculate_loc_Y(X_mean[,1], X_mean[,2], D_range)

new_data <- cbind(X_mean, D_range)
predicted_y_CATE <- predict(dbfit, X_test= cbind(X_mean[,1], X_mean[,2]), Z_test = D_range,pi_test = X_mean[,3])
all_estimated_theoretical_y <- predicted_y_CATE$y_hat[,sample_inds]
estimated_theoretical_y <- rowMeans(all_estimated_theoretical_y)

theoretical_CATE_RMSE[j,i] <- sqrt(mean((loc_Y_D - estimated_theoretical_y)^2))
theoretical_CATE_MAE[j,i] <- mean(abs(loc_Y_D - estimated_theoretical_y))
theoretical_CATE_MAPE[j,i] <- mean(abs(loc_Y_D - estimated_theoretical_y)/abs(loc_Y_D))

# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
theoretical_CATE_lower_bounds <- apply(all_estimated_theoretical_y, 1, function(x) quantile(x, probs = 0.025))
theoretical_CATE_upper_bounds <- apply(all_estimated_theoretical_y, 1, function(x) quantile(x, probs = 0.975))

theoretical_CATE_len[j,i] <- mean(abs(theoretical_CATE_upper_bounds-theoretical_CATE_lower_bounds))

cover_vector <- (loc_Y_D >= theoretical_CATE_lower_bounds) & (loc_Y_D <= theoretical_CATE_upper_bounds)

# Calculate the number of times tau is within the CI
number_covered <- sum(cover_vector)

theoretical_CATE_cover[j,i] <- number_covered / length(loc_Y_D)

# Define the y-range to cover all values
y_min <- min(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
y_max <- max(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))

png("BCF_CATE.png", width = 800, height = 600)
plot(D_range, estimated_theoretical_y, 
     main = "Estimated CATE (BCF)", 
     xlab = "D", 
     ylab = "Y",
     pch = 16,  # solid circles
     col = "blue",
     ylim = c(y_min, y_max))

# Add dotted lines for the lower and upper bounds
lines(D_range, theoretical_CATE_lower_bounds, type = "l", col = "blue", lwd = 1, lty = 2)
lines(D_range, theoretical_CATE_upper_bounds, type = "l", col = "blue", lwd = 1, lty = 2)

# Add the line plot of D_range and loc_Y_D
lines(D_range, loc_Y_D, type = "l", col = "red", lwd = 2)

# Add a legend to distinguish between the scatter plot and the line plot
legend("topright", legend = c("Estimated", "Actual", "Lower Bound", "Upper Bound"), 
       col = c("blue", "red", "blue", "blue"), pch = c(16, NA, NA, NA), 
       lty = c(NA, 1, 2, 2), lwd = c(NA, 2, 1, 1))
dev.off()

}

}