coef_valid <- function(
    X,
    coef_R,
    coef_Cpp,
    tolerance = 1e-6) {
  
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

sens_plot <- function(
    X,
    Y,
    eta = seq(0.1, 0.9, 0.1),
    K = 1) {
  
  data <- matrix(nrow = 0, ncol = 3)
  
  for (i in 1:length(eta)) {
    spls_mod <- spls_cpp(
      X,
      Y,
      K = K,
      eta = eta[i],
      kappa = 0.5,
      select = "pls2",
      fit = "widekernelpls",
      scale_x = TRUE,
      scale_y = FALSE,
      eps = 1e-4,
      maxstep = 100,
      trace = FALSE)
    
    coef <- spls_mod$betahat
    coef_idx <- which(coef != 0)
    coef <- coef[coef_idx]
    print(coef)
    ids <- colnames(X)[coef_idx]
    data <- rbind(
      data,
      matrix(cbind(ids, abs(coef), rep(eta[i], length(coef))), ncol = 3)
    )
  }
  
  data <- data %>%
    as.data.frame() %>%
    rename("SNP_ID" = "V1",
           "Abs_Coefficient" = "V2",
           "eta" = "V3") %>%
    mutate(Abs_Coefficient = as.numeric(Abs_Coefficient),
           SNP_ID = as.character(SNP_ID),
           eta = as.numeric(eta))
  
  ggplot(data, aes(x = SNP_ID, y = Abs_Coefficient, color = factor(eta))) +
    geom_point(size = 3, alpha = 0.8) +
    coord_cartesian(ylim = c(
      min(floor(data$Abs_Coefficient * 100) / 100),
      max(ceiling(data$Abs_Coefficient * 100) / 100))) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = expression("Coefficients for Selected SNPs by " * eta),
         x = "SNP ID",
         y = expression("|"~hat(beta)~"|"),
         color = expression(eta)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}