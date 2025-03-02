library(Rcpp)
# Load the Rcpp-based PLS function
sourceCpp("plsr_rcpp.cpp")
source("plsr_rcpp.R")

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

