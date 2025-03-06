#ifndef SPLS_DV_H
#define SPLS_DV_H

#include <RcppArmadillo.h>

arma::vec ust(const arma::vec& b, double eta);
arma::vec spls_dv(const arma::mat& Z, double eta, double kappa, double eps, int maxstep);

#endif // SPLS_DV_H