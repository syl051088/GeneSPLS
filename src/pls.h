#ifndef PLS_H
#define PLS_H

#include <RcppArmadillo.h>
  
Rcpp::List widekernelpls_fit(arma::mat X,
                               arma::mat Y,
                               int ncomp = 5,
                               bool center = true,
                               double tol = 1.5e-8,
                               int maxit = 100);
 
 
#endif // PLS_H   
   
   