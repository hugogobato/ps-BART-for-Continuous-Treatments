# Number of samples per dataset
num_samples <- 100

set.seed(42)  # For reproducibility

for (i in 1:100) {
  X1 <- rexp(num_samples, rate=1)
  X2 <- rexp(num_samples, rate=1)

  D <- rexp(num_samples, rate=1 / (X1 + X2))

  Y <- rnorm(num_samples, mean = D + (X1 + X2) * exp(-D * (X1 + X2)), sd = 1)
  E_Y <- D + (X1 + X2) * exp(-D * (X1 + X2))

  # Create a data frame for the dataset
  data <- data.frame(
    X1 = X1,
    X2 = X2,
    D = D,
    Y = Y,
    E_Y = E_Y
  )

  # Save the dataset as a CSV file
  file_name <- paste0("dataset_", i, ".csv")
  write.csv(data, file = file_name, row.names = FALSE)
}

