#ifndef SPLS_H
#define SPLS_H

#include <RcppArmadillo.h>

Rcpp::List spls_cpp(arma::mat x, arma::mat y, int K, double eta, double kappa,
                    std::string select, std::string fit,
                    bool scale_x, bool scale_y, 
                    double eps, int maxstep, bool trace);

#endif // SPLS_H