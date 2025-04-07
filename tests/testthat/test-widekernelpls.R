library(testthat)
library(pls)
library(GeneSPLS)

# Test function for comparing projection and coefficients with multivariate Y
test_that("Projection and coefficients match between R and C++ for multivariate Y", {
  # Set seed for reproducibility
  set.seed(815)
  
  # Create simulated data with p > n and multivariate Y
  n <- 20    # Number of observations
  p <- 50    # Number of predictors (p > n)
  m <- 3     # Number of response variables
  
  # Create predictor matrix X
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Create response matrix Y
  Y <- matrix(0, nrow = n, ncol = m)
  Y[,1] <- X[,1:5] %*% c(1, 0.5, -0.5, -1, 0.8) + rnorm(n, 0, 0.2)
  Y[,2] <- X[,3:7] %*% c(0.7, 1.2, -1, 0.5, -0.3) + rnorm(n, 0, 0.2)
  Y[,3] <- X[,6:10] %*% c(-0.8, 1, 0.6, -0.2, 0.9) + rnorm(n, 0, 0.2)
  
  # Create data frame for R implementation
  sim_data <- data.frame(Y1 = Y[,1], Y2 = Y[,2], Y3 = Y[,3], X)
  
  # Number of components
  ncomp <- 5
  
  # Run R implementation
  pls_r <- pls::plsr(cbind(Y1, Y2, Y3) ~ ., 
                data = sim_data, 
                ncomp = ncomp, 
                method = "widekernelpls", 
                validation = "none",
                center = TRUE)
  
  # Run Rcpp implementation
  pls_cpp <- widekernelpls_fit(X, Y, ncomp, center = TRUE)
  
  # PART 1: COMPARE PROJECTION MATRICES
  # Extract projection matrices
  proj_r <- pls_r$projection
  proj_cpp <- pls_cpp$projection
  
  # Test 1: Projection matrices should have the same dimensions
  expect_equal(dim(proj_r), dim(proj_cpp), 
               info = "Projection matrices have different dimensions")
  
  # Test 2: Correlation between projection matrices should be very high
  # Use absolute correlation since signs might be flipped but still valid
  proj_cor <- abs(cor(as.vector(proj_r), as.vector(proj_cpp)))
  expect_gte(proj_cor, 0.99)
  
  # PART 2: COMPARE COEFFICIENTS
  # Extract coefficients for each component number (adjusting for format differences)
  for (nc in 1:ncomp) {
    # R plsr returns coefficients as [p, m, 1] for each ncomp
    # C++ returns coefficients as [p, m, ncomp] with a slice for each ncomp
    coef_r <- coef(pls_r, ncomp = nc)[,,1]
    coef_cpp <- pls_cpp$coefficients[,,nc]
    
    # Test 3: Coefficients should have the same dimensions
    expect_equal(dim(coef_r), dim(coef_cpp), 
                 info = paste("Component", nc, "coefficients have different dimensions"))
    
    # Test 4: For each response variable, coefficients should be highly correlated
    for (i in 1:m) {
      coef_r_i <- coef_r[,i]
      coef_cpp_i <- coef_cpp[,i]
      
      # Correlation should be high (use absolute value to account for sign flips)
      coef_cor <- abs(cor(coef_r_i, coef_cpp_i))
      
      expect_gte(coef_cor, 0.99)
    }
  }
})