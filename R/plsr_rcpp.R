# Load Rcpp
library(Rcpp)

# Compile and load the C++ implementation
sourceCpp("plsr_rcpp.cpp")

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
