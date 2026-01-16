set.seed(123)  # For reproducibility

# Function to generate confounders
generate_confounders <- function(N) {
  C1 <- rnorm(N)
  C2 <- rnorm(N)
  C3 <- rnorm(N)
  C4 <- rnorm(N)
  C5 <- runif(N, -2, 2)
  C6 <- runif(N, -3, 3)
  data.frame(C1, C2, C3, C4, C5, C6)
}

# Function to generate W based on the specified model
generate_W <- function(C, spec) {
  if (spec == 1) {
    W <- 9 * (0.8 + 0.1 * C$C1 + 0.1 * C$C2 - 0.1 * C$C3 + 0.2 * C$C4 + 0.1 * C$C5) - 3 + rnorm(nrow(C), 0, 5)
  } else if (spec == 2) {
    W <- 15 * (0.8 + 0.1 * C$C1 + 0.1 * C$C2 - 0.1 * C$C3 + 0.2 * C$C4 + 0.1 * C$C5) + 2 + 2 * rt(nrow(C), df = 4)
  } else if (spec == 3) {
    W <- 15 * (0.8 + 0.1 * C$C1 + 0.1 * C$C2 - 0.1 * C$C3 + 0.2 * C$C4 + 0.1 * C$C5) + 3 * sqrt(2) * C$C3^2 + rt(nrow(C), df = 4)
  }
  return(W)
}

# Function to generate Y based on W and C
generate_Y_and_mu <- function(W, C) {
  mu_WC <- -10 - (2 * C$C1 + 2 * C$C2 + 3 * C$C3 - 1 * C$C4 + 2 * C$C5) - 
    W * (0.1 - 0.1 * C$C1 + 0.1 * C$C4 + 0.1 * C$C5 + 0.1 * C$C3^2) + 0.13 * W^3
  Y <- rnorm(length(W), mean = mu_WC, sd = sqrt(10))
  return(list(Y = Y, mu_WC = mu_WC))
}

# Parameters
sample_sizes <- c(100, 250, 500)
num_datasets <- 100

# Initialize the list to store datasets
data_list <- list()

# Generate datasets
for (N in sample_sizes) {
  for (dataset_num in 1:num_datasets) {
    for (spec in 1:3) {
      # Generate confounders
      confounders <- generate_confounders(N)
      
      # Generate W
      W <- generate_W(confounders, spec)
      
       # Generate Y and mu_WC
      result <- generate_Y_and_mu(W, confounders)
      Y <- result$Y
      mu_WC <- result$mu_WC
      
      # Combine into a single data frame
      data <- cbind(confounders, W, Y, mu_WC)
      
      # Save or use the dataset (e.g., save to file, store in list, etc.)
      # For example, saving to a list (assuming data_list is a predefined empty list):
       data_list[[paste0("N", N, "_dataset", dataset_num, "_spec", spec)]] <- data
      
      # Uncomment the following line if you want to save each dataset to a CSV file:
       write.csv(data, file = paste0("data_N", N, "_dataset", dataset_num, "_spec", spec, ".csv"), row.names = FALSE)
    }
  }
}

# Print message to indicate completion
cat("Datasets generated and saved successfully.\n")
