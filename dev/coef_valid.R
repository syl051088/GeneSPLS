coef_valid <- function(X, coef_R, coef_Cpp, tolerance = 1e-6) {
  
  # Identify non-zero coefficients using linear indices
  nonzero_idx_R <- which(coef_R != 0)
  nonzero_idx_Cpp <- which(coef_Cpp != 0)
  both <- intersect(nonzero_idx_R, nonzero_idx_Cpp)
  SNPs_Both <- colnames(X)[both]
  
  cat("Number of non-zero coefficients (spls):", length(nonzero_idx_R), "\n")
  cat("Number of non-zero coefficients (GeneSPLS):", length(nonzero_idx_Cpp), "\n")
  
  # Check if the selected SNPs (non-zero coefficient indices) match exactly
  if(identical(nonzero_idx_R, nonzero_idx_Cpp)) {
    
    cat("Selected SNPs are identical between the two models.\n")
    cat("SNPs selected:", SNPs_Both, "\n")
    
    max_diff <- max(abs(coef_R[nonzero_idx_R] - coef_Cpp[nonzero_idx_R]))
    
    if(max_diff < tolerance) {
      
      cat("Both sets of coefficients are equal within a tolerance of", tolerance, "\n")
      
    } else {
      
      cat("Warning: Maximum difference in corresponding coefficients is", max_diff,
          "which exceeds the tolerance of", tolerance, "\n")
    }
    
  } else {
    
    cat("Warning: selected SNPs differ between the two models.\n")
    
    # Report differences in indices
    diff_in_R <- setdiff(nonzero_idx_R, nonzero_idx_Cpp)
    diff_in_Cpp <- setdiff(nonzero_idx_Cpp, nonzero_idx_R)
    SNPs_R <- colnames(X)[diff_in_R]
    SNPs_Cpp <- colnames(X)[diff_in_Cpp]
    
    if(length(diff_in_R) > 0) {
      
      cat("SNPs selected by both:", SNPs_Both, "\n")
      cat("SNPs selected by spls but not by GeneSPLS:", SNPs_R, "\n")
      
    } 
    
    if(length(diff_in_Cpp) > 0) {
      
      cat("SNPs selected by both:", SNPs_Both, "\n")
      cat("SNPs selected by GeneSPLS but not by spls:", SNPs_Cpp, "\n")
      
    }
  }
}