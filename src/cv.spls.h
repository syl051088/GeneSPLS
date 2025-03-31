#ifndef CV_SPLS_H
#define CV_SPLS_H

#include <RcppArmadillo.h>

// Declaration of the cv_spls_cpp function
Rcpp::List cv_spls_cpp(arma::mat x,
                       arma::mat y,
                       int fold,
                       Rcpp::NumericVector eta,
                       Rcpp::IntegerVector K,
                       double kappa = 0.5,
                       std::string select = "pls2",
                       bool scale_x = true,
                       bool scale_y = false,
                       double eps = 1e-4,
                       int maxstep = 100,
                       bool trace = false);

#endif // CV_SPLS_H
