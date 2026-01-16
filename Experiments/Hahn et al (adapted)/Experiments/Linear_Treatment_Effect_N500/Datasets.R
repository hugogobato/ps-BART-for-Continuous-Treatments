library(dplyr)
library(purrr)
library(readr)

# Define the provided functions
g <- function(x) {
  case_when(
    x == 1 ~ 2,
    x == 2 ~ -1,
    x == 3 ~ -4
  )
}

generate_x_values <- function(n) {
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n, 0, 1)
  x3 <- rnorm(n, 0, 1)
  x4 <- sample(1:3, n, replace = TRUE)
  x5 <- sample(0:1, n, replace = TRUE)
  list(x1, x2, x3, x4, x5)
}

generate_D_values <- function(n, x1, x5, x2, x3) {
  alpha <- 2
  D <- rexp(n, rate = alpha / (abs(x1) + x5 + abs(x2 - x3)))
  return(D)
}

simulate <- function(n, tau_func, mu_func) {
  values <- generate_x_values(n)
  names(values) <- c("x1", "x2", "x3", "x4", "x5")
  x4_mapped <- map_dbl(values$x4, g)
  
  mu <- mu_func(x4_mapped, values$x1, values$x3)
  tau <- tau_func(values$x2, values$x5)
  treatment <- generate_D_values(n, values$x1, values$x5, values$x2, values$x3)
  
  outcome <- mu + treatment * tau 
  
  data.frame(x1 = values$x1, x2 = values$x2, x3 = values$x3, x4 = values$x4, x5 = values$x5, 
             D = treatment, y = outcome, tau = tau, mu = mu)
}

simulate_heterogeneous_nonlinear <- function(n) {
  tau_func <- function(x2, x5) { 1 + 2 * x2 * x5 }
  mu_func <- function(x4, x1, x3) { -6 + x4 + 6 * abs(x3 - 1) }
  simulate(n, tau_func, mu_func)
}

# Define the number of iterations
num_iterations <- 100
n <- 500  # Define the number of samples per dataset

# Generate and save datasets
for (i in 1:num_iterations) {
  dataset <- simulate_heterogeneous_nonlinear(n)
  write_csv(dataset, paste0("dataset_", i, ".csv"))
}
