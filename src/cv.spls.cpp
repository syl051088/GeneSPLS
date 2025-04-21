#include "cv.spls.h"
#include "spls.h"       // for spls_cpp()
#include "correctp.h"   // for correctp()
#include <algorithm>
#include <vector>
#include <random>
#include <ctime>

// [[Rcpp::depends(RcppArmadillo)]]

// Utility function to generate CV folds (returns a vector of index vectors)
std::vector< std::vector<int> > generate_folds(int n, int fold) {
  std::vector<int> indices(n);
  for (int i = 0; i < n; i++) {
    indices[i] = i;
  }
  
  // Use a random device and engine for shuffling
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);
  
  std::vector< std::vector<int> > folds(fold);
  for (int i = 0; i < n; i++) {
    folds[i % fold].push_back(indices[i]);
  }
  return folds;
}


//' @name cv_spls_cpp
//' @title Cross-Validation for SPLS
//'
//' @description This function performs cross-validation for the number of
//' latent components and sparsity parameters
//'
//' @param x Predictor matrix (n x p)
//' @param y Response matrix (n x q)
//' @param fold Number of folds for cross-validation (default: 5)
//' @param eta Numeric vector of sparsity parameters (default: (0.95, 0.96, 0.97, 0.98, 0.99))
//' @param K Integer vector of candidate numbers of components (default: (1, 2, 3, 4, 5))
//' @param kappa Sparsity parameter (default: 0.5)
//' @param select Selection method (default: "pls2")
//' @param scale_x Logical, whether to scale x (default: TRUE)
//' @param scale_y Logical, whether to scale y (default: FALSE)
//' @param eps Convergence criterion (default: 1e-4)
//' @param maxstep Maximum number of iterations (default: 100)
//' @param trace Logical, whether to print progress (default: FALSE)
//'
//' @return A list containing:
//' \item{mspemat}{Matrix of mean squared prediction errors (MSPE)}
//' \item{eta.opt}{Optimal eta value}
//' \item{K.opt}{Optimal number of latent components}
//' 
//' @useDynLib GeneSPLS, .registration = TRUE
//' @import Rcpp
//' 
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List cv_spls_cpp(arma::mat x,
                       arma::mat y,
                       int fold,
                       Rcpp::NumericVector eta,
                       Rcpp::IntegerVector K,
                       double kappa,
                       std::string select,
                       bool scale_x,
                       bool scale_y,
                       double eps,
                       int maxstep,
                       bool trace) {
  
  // Basic input checks
  if(x.n_rows != y.n_rows) {
    Rcpp::stop("Number of rows in x and y must be equal.");
  }
  
  int n = x.n_rows; // number of observations
  int p = x.n_cols;
  int q = y.n_cols;
  
  // Convert candidate vectors to std::vector for convenience
  std::vector<double> eta_candidates = Rcpp::as< std::vector<double> >(eta);
  std::vector<int> K_candidates = Rcpp::as< std::vector<int> >(K);
  
  int n_eta = eta_candidates.size();
  int n_K = K_candidates.size();
  
  // Create an MSPE matrix (eta x K) to store the averaged errors.
  arma::mat mspemat(n_eta, n_K, arma::fill::zeros);
  
  // Generate CV folds (each fold is a vector of row indices)
  std::vector< std::vector<int> > fold_indices = generate_folds(n, fold);
  
  // Loop over candidate eta values
  for (int i = 0; i < n_eta; i++) {
    double current_eta = eta_candidates[i];
    Rcpp::Rcout << "eta = " << current_eta << "\n";
    
    // Matrix to store MSPE for each fold (rows = folds, columns = candidate K values)
    arma::mat mspemati(fold, n_K, arma::fill::zeros);
    
    // Loop over CV folds
    for (int j = 0; j < fold; j++) {
      // Determine test indices for this fold
      std::vector<int> test_idx = fold_indices[j];
      
      // Determine training indices (all indices not in test_idx)
      std::vector<int> train_idx;
      train_idx.reserve(n - test_idx.size());
      // Create a boolean mask for test indices
      std::vector<bool> is_test(n, false);
      for (int idx : test_idx) {
        is_test[idx] = true;
      }
      for (int idx = 0; idx < n; idx++) {
        if (!is_test[idx]) {
          train_idx.push_back(idx);
        }
      }
      
      // Create training and test matrices
      arma::mat x_train(train_idx.size(), p);
      arma::mat y_train(train_idx.size(), q);
      arma::mat x_test(test_idx.size(), p);
      arma::mat y_test(test_idx.size(), q);
      
      for (size_t r = 0; r < train_idx.size(); r++) {
        x_train.row(r) = x.row(train_idx[r]);
        y_train.row(r) = y.row(train_idx[r]);
      }
      for (size_t r = 0; r < test_idx.size(); r++) {
        x_test.row(r) = x.row(test_idx[r]);
        y_test.row(r) = y.row(test_idx[r]);
      }
      
      // Call correctp() to validate and possibly adjust parameters
      Rcpp::List cp = correctp(x_train, y_train, current_eta, K_candidates[0], kappa, select, "widekernelpls");
      // Extract corrected parameters (if needed, though in our case they may remain unchanged)
      double corrected_eta = Rcpp::as<double>(cp["eta"]);
      // int corrected_K_dummy = Rcpp::as<int>(cp["K"]);  // not used directly; we use candidate K values
      kappa = Rcpp::as<double>(cp["kappa"]);
      select = Rcpp::as<std::string>(cp["select"]);
      std::string fit = "widekernelpls"; // force to "widekernelpls"
      
      // For the CV step, we call spls_cpp() with K = max(candidate K)
      int K_max = *std::max_element(K_candidates.begin(), K_candidates.end());
      Rcpp::List spls_obj = spls_cpp(x_train, y_train, K_max, corrected_eta, kappa, select,
                                     fit, scale_x, scale_y, eps, maxstep, trace);
      
      // Extract required components from spls_obj
      Rcpp::List betamat = spls_obj["betamat"]; // list of coefficient matrices, index 0 corresponds to 1 component, etc.
      arma::rowvec meanx = spls_obj["meanx"];
      arma::rowvec normx = spls_obj["normx"];
      arma::rowvec mu = spls_obj["mu"];
      
      // Scale x_test using training mean and norm (replicating R's scale() function)
      arma::mat newx = x_test;
      for (unsigned int col = 0; col < newx.n_cols; col++) {
        newx.col(col) = (newx.col(col) - meanx(col)) / normx(col);
      }
      
      // For each candidate K, compute predictions and MSPE on the test fold
      for (int k_ind = 0; k_ind < n_K; k_ind++) {
        int current_K = K_candidates[k_ind];
        // In spls_cpp, betamat is a list of length K_max, where betamat[i] corresponds to i+1 latent components.
        // Use the coefficient matrix for current_K (i.e. index current_K - 1)
        arma::mat beta = Rcpp::as<arma::mat>(betamat[current_K - 1]);
        
        // Compute prediction: pred = newx * beta + (vector of ones) * mu
        arma::vec ones = arma::ones<arma::vec>(newx.n_rows);
        arma::mat pred = newx * beta;
        // Add intercept mu (assumed to be a row vector; add as each row)
        for (unsigned int r = 0; r < pred.n_rows; r++) {
          pred.row(r) = pred.row(r) + mu;
        }
        
        // Compute mean squared error: average over observations and responses
        arma::mat diff = y_test - pred;
        double mse = arma::accu(arma::square(diff)) / (diff.n_rows * diff.n_cols);
        mspemati(j, k_ind) = mse;
      } // end candidate K loop
      
    } // end fold loop
    
    // Average MSPE over folds for current eta candidate (for each candidate K)
    for (int k_ind = 0; k_ind < n_K; k_ind++) {
      double mean_mspe = arma::mean(mspemati.col(k_ind));
      mspemat(i, k_ind) = mean_mspe;
    }
    
  } // end eta loop
  
  // Determine optimal parameters: 
  // minpmse is the minimum MSPE over all candidate combinations
  double minpmse = mspemat.min();
  
  // Compute row and column minima
  arma::colvec rowMin = arma::min(mspemat, 1);  // minimum of each row (for each eta)
  arma::rowvec colMin = arma::min(mspemat, 0);    // minimum of each column (for each candidate K)
  
  // Optimal K: among candidate K where the column minimum equals minpmse, choose the smallest K.
  int K_opt = 0;
  bool foundK = false;
  for (int k_ind = 0; k_ind < n_K; k_ind++) {
    if (std::abs(colMin(k_ind) - minpmse) < 1e-8) { // accounting for floating point precision
      if (!foundK || K_candidates[k_ind] < K_opt) {
        K_opt = K_candidates[k_ind];
        foundK = true;
      }
    }
  }
  
  // Optimal eta: among candidate eta where the row minimum equals minpmse, choose the maximum eta.
  double eta_opt = 0;
  bool foundEta = false;
  for (int i = 0; i < n_eta; i++) {
    if (std::abs(rowMin(i) - minpmse) < 1e-8) {
      if (!foundEta || eta_candidates[i] > eta_opt) {
        eta_opt = eta_candidates[i];
        foundEta = true;
      }
    }
  }
  
  Rcpp::Rcout << "\nOptimal parameters: eta = " << eta_opt << ", K = " << K_opt << "\n";
  
  // Return the results in a list similar to the R version
  return Rcpp::List::create(Rcpp::Named("mspemat") = mspemat,
                            Rcpp::Named("eta.opt") = eta_opt,
                            Rcpp::Named("K.opt") = K_opt);
}
