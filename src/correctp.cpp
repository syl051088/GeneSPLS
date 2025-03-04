#include <RcppArmadillo.h>
#include <algorithm>  // Required for std::find

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List correctp(arma::mat X, arma::mat Y, double eta, int K, double kappa,
                    std::string select, std::string fit) {

  // Validate eta
  if (eta < 0 || eta >= 1) {
    stop("eta should be strictly between 0 and 1.");
  }

  // Validate K
  if (K > (int)X.n_cols) {
    stop("K cannot exceed the number of predictors.");
  }
  if (K >= (int)X.n_rows) {
    stop("K cannot exceed the sample size.");
  }
  if (K <= 0) {
    stop("K should be a positive integer.");
  }

  // Validate kappa
  if (kappa > 0.5 || kappa < 0) {
    Rcout << "Warning: kappa should be between 0 and 0.5. Setting kappa=0.5.\n";
    kappa = 0.5;
  }

  // Validate selection method
  if (select != "pls2" && select != "simpls") {
    Rcout << "Invalid PLS algorithm for variable selection. Using 'pls2'.\n";
    select = "pls2";
  }

  // Validate PLS fitting method
  std::vector<std::string> valid_fits = {"simpls", "kernelpls", "widekernelpls", "oscorespls"};
  if (std::find(valid_fits.begin(), valid_fits.end(), fit) == valid_fits.end()) {
    Rcout << "Invalid PLS algorithm for model fitting. Using 'simpls'.\n";
    fit = "simpls";
  }

  return List::create(
    Named("K") = K,
    Named("eta") = eta,
    Named("kappa") = kappa,
    Named("select") = select,
    Named("fit") = fit
  );
}
