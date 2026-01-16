
number_of_points <- 500
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

i=1

file_name <- paste0("dataset_", i, ".csv")
data<-read.csv(file_name)
Y <- data$y
X <- data[, !(names(data) %in% c("y", "tau", "mu", "D"))]
D <- data$D

max_D <- max(D)
min_D <- min(D)

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

real_ATE <-  calculate_loc_Y(X_mean[,1], X_mean[,2],X_mean[,3], X_mean[,4],X_mean[,5], D)

# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
empirical_ATE_lower_bounds <- apply(bart_model$yhat.test, 2, function(x) quantile(x, probs = 0.025))
empirical_ATE_upper_bounds <- apply(bart_model$yhat.test, 2, function(x) quantile(x, probs = 0.975))

D_range <- seq(min_D, max_D, length.out = number_of_points)

loc_Y_D <- calculate_loc_Y(X_mean[,1], X_mean[,2],X_mean[,3], X_mean[,4],X_mean[,5], D_range)
# Define the y-range to cover all values
y_min <- min(c(predicted_y, empirical_ATE_lower_bounds, empirical_ATE_upper_bounds, loc_Y_D))
y_max <- max(c(predicted_y, empirical_ATE_lower_bounds, empirical_ATE_upper_bounds, loc_Y_D))
png("ps-BART_observed_D.png", width = 800, height = 600)
plot(D, predicted_y, 
     main = "Observed D (ps-BART)", 
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
all_estimated_theoretical_y <- predict(bart_model, newdata = new_data,)
estimated_theoretical_y <- colMeans(all_estimated_theoretical_y)
# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
theoretical_CATE_lower_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.025))
theoretical_CATE_upper_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.975))
# Define the y-range to cover all values
y_min <- min(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
y_max <- max(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
png("ps-BART_ATE.png", width = 800, height = 600)
plot(D_range, estimated_theoretical_y, 
     main = "Estimated ATE (ps-BART)", 
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
j=1
column_means[1] <- X$x1[j]
column_means[2] <- X$x2[j]
column_means[3] <- X$x3[j]
column_means[4] <- X$x4[j]
column_means[5] <- X$x5[j]
column_means[6] <- X$pi_hats[j]

X_mean <- matrix(rep(column_means, each = number_of_points), nrow = number_of_points, byrow = FALSE)
loc_Y_D <- calculate_loc_Y(X_mean[,1], X_mean[,2],X_mean[,3], X_mean[,4],X_mean[,5], D_range)

new_data <- cbind(X_mean, D_range)
all_estimated_theoretical_y <- predict(bart_model, newdata = new_data,)
estimated_theoretical_y <- colMeans(all_estimated_theoretical_y)
# Calculate the 2.5th and 97.5th percentiles for each datapoint across simulations
theoretical_CATE_lower_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.025))
theoretical_CATE_upper_bounds <- apply(all_estimated_theoretical_y, 2, function(x) quantile(x, probs = 0.975))
# Define the y-range to cover all values
y_min <- min(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
y_max <- max(c(estimated_theoretical_y, theoretical_CATE_lower_bounds, theoretical_CATE_upper_bounds, loc_Y_D))
png("ps-BART_CATE.png", width = 800, height = 600)
plot(D_range, estimated_theoretical_y, 
     main = "Estimated CATE (ps-BART)", 
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