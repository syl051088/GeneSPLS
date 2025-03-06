#include "correctp.h"
#include <algorithm>  // Required for std::find
// [[Rcpp::depends(RcppArmadillo)]]

//' @name correctp
//' @title Parameter validation for SPLS
//' 
//' @description Validates and corrects parameters for the SPLS algorithm
//' 
//' @param X Predictor matrix
//' @param Y Response matrix
//' @param eta Sparsity parameter between 0 and 1
//' @param K Number of components
//' @param kappa Parameter between 0 and 0.5
//' @param select Selection method ("pls2" or "simpls")
//' @param fit Fitting method
//' @return List with validated parameters
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List correctp(arma::mat X, arma::mat Y, double eta, int K, double kappa,
                   std::string select, std::string fit) {
 
 // Validate eta
 if (eta < 0 || eta >= 1) {
   Rcpp::stop("eta should be strictly between 0 and 1.");
 }
 
 // Validate K
 if (K > (int)X.n_cols) {
   Rcpp::stop("K cannot exceed the number of predictors.");
 }
 if (K >= (int)X.n_rows) {
   Rcpp::stop("K cannot exceed the sample size.");
 }
 if (K <= 0) {
   Rcpp::stop("K should be a positive integer.");
 }
 
 // Validate kappa
 if (kappa > 0.5 || kappa < 0) {
   Rcpp::Rcout << "Warning: kappa should be between 0 and 0.5. Setting kappa=0.5.\n";
   kappa = 0.5;
 }
 
 // Validate selection method
 if (select != "pls2" && select != "simpls") {
   Rcpp::Rcout << "Invalid PLS algorithm for variable selection. Using 'pls2'.\n";
   select = "pls2";
 }
 
 // Validate PLS fitting method - always use widekernelpls for this implementation
 std::vector<std::string> valid_fits = {"simpls", "kernelpls", "widekernelpls", "oscorespls"};
 if (std::find(valid_fits.begin(), valid_fits.end(), fit) == valid_fits.end()) {
   Rcpp::Rcout << "Invalid PLS algorithm for model fitting. Using 'widekernelpls'.\n";
 }
 
 // Always set to widekernelpls
 fit = "widekernelpls";
 
 return Rcpp::List::create(
   Rcpp::Named("K") = K,
   Rcpp::Named("eta") = eta,
   Rcpp::Named("kappa") = kappa,
   Rcpp::Named("select") = select,
   Rcpp::Named("fit") = fit
 );
}