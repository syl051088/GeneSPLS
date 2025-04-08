#ifndef SPLS_H
#define SPLS_H

#include <RcppArmadillo.h>

Rcpp::List spls_cpp(arma::mat x,
                    arma::mat y,
                    int K = 1,
                    double eta = 0.99,
                    double kappa = 0.5,
                    std::string select = "pls2",
                    std::string fit = "widekernelpls",
                    bool scale_x = true,
                    bool scale_y = false, 
                    double eps = 1e-4,
                    int maxstep = 100,
                    bool trace = false);

#endif // SPLS_H