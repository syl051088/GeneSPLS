#ifndef CV_SPLS_H
#define CV_SPLS_H

#include <RcppArmadillo.h>

// Declaration of the cv_spls_cpp function
Rcpp::List cv_spls_cpp(arma::mat x,
                       arma::mat y,
                       int fold = 5,
                       Rcpp::NumericVector eta = Rcpp::NumericVector::create(0.95, 0.96, 0.97, 0.98, 0.99),
                       Rcpp::IntegerVector K = Rcpp::IntegerVector::create(1, 2, 3, 4, 5),
                       double kappa = 0.5,
                       std::string select = "pls2",
                       bool scale_x = true,
                       bool scale_y = false,
                       double eps = 1e-4,
                       int maxstep = 100,
                       bool trace = false);

#endif // CV_SPLS_H
