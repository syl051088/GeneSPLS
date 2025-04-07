#ifndef PLS_H
#define PLS_H

#include <RcppArmadillo.h>
  
Rcpp::List widekernelpls_fit(arma::mat X,
                               arma::mat Y,
                               int ncomp,
                               bool center,
                               double tol,
                               int maxit);
 
 
#endif // PLS_H   
   
   