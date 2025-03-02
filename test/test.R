library(Rcpp)
library(GeneSPLS)

# Define the wrapper function
plsr_rcpp <- function(Y, X, ncomp = 5) {
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)
  
  # Ensure Y is correctly formatted as a matrix
  if (ncol(Y_mat) == 1) {
    Y_mat <- matrix(Y_mat, ncol = 1)
  }
  
  # Call the C++ function
  model <- widekernelpls_rcpp(X_mat, Y_mat, ncomp)
  
  # Format results as a list
  result <- list(
    coefficients = model$coefficients,
    scores = model$scores,
    loadings = model$loadings,
    projection = model$projection,
    fitted.values = model$fitted_values,
    residuals = model$residuals
  )
  
  class(result) <- "plsr_rcpp"
  return(result)
}

# Summary function for easy printing
summary.plsr_rcpp <- function(object) {
  cat("PLSR Model Summary:\n")
  print(object)
}

# Simulated Data
set.seed(42)
n <- 20
p <- 100
X <- matrix(rnorm(n * p), n, p)

# Generate proper Y matrix
beta_true <- rep(0, p)
beta_true[sample(1:p, 5)] <- runif(5, -2, 2)
Y <- as.matrix(X %*% beta_true + rnorm(n, sd = 0.5))  # Ensure Y is a matrix

# Check dimensions before running
print(dim(Y))  # Should return [20, 1]

# Fit PLSR Model
pls_model_star <- plsr_rcpp(Y, X, ncomp = 5)

# Print Summary
summary(pls_model_star)

