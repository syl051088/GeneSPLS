#ifndef CORRECTP_H
#define CORRECTP_H

#include <RcppArmadillo.h>

Rcpp::List correctp(arma::mat X, arma::mat Y, double eta, int K, double kappa,
                     std::string select, std::string fit);

#endif // CORRECTP_H