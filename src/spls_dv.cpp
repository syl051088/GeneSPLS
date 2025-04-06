#include "spls_dv.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @name ust
//' @title Univariate Soft Thresholding (UST) Function
//' 
//' @description Performs soft thresholding on a vector
//' 
//' @param b Vector to be thresholded
//' @param eta Threshold parameter between 0 and 1
//' @return Soft thresholded vector
//' 
arma::vec ust(const arma::vec& b, double eta) {
 arma::vec b_ust = arma::zeros<arma::vec>(b.n_elem);
 if (eta < 1) {
   arma::vec valb = arma::abs(b) - eta * arma::max(arma::abs(b));
   arma::uvec idx = arma::find(valb >= 0);
   b_ust.elem(idx) = valb.elem(idx) % arma::sign(b.elem(idx));
 }
 return b_ust;
}

//' @name spls_dv
//' @title Sparse PLS Direction Vector Estimation
//' 
//' @description Calculates the direction vector for Sparse PLS
//' 
//' @param Z Cross-product matrix (p x q)
//' @param eta Sparsity parameter between 0 and 1
//' @param kappa Parameter between 0 and 0.5
//' @param eps Convergence criterion
//' @param maxstep Maximum number of iterations
//' @return Direction vector
//' 
arma::vec spls_dv(const arma::mat& Z, double eta, double kappa, double eps, int maxstep) {
 int p = Z.n_rows;
 int q = Z.n_cols;
 double Znorm1 = arma::median(arma::abs(arma::vectorise(Z)));
 arma::mat Z_scaled = Z / Znorm1;
 
 // Case: Univariate response
 if(q == 1) {
   return ust(Z_scaled.col(0), eta);
 }
 
 // Case: Multivariate response
 arma::mat M = Z_scaled * Z_scaled.t();  // (p x p)
 
 // Initialize c as in the original R code: (p x 1) vector with constant value 10
 arma::vec c = 10 * arma::ones<arma::vec>(p);
 arma::vec c_old = c;
 double dis = 10.0;
 int i = 1;
 
 if(kappa == 0.5) {  // Use SVD if kappa == 0.5
   while(dis > eps && i <= maxstep) {
     // Compute M * c as a vector and force it to be a (p x 1) matrix
     arma::vec temp = M * c;
     arma::mat M_c(temp);
     
     // Perform SVD on the (p x 1) matrix
     arma::mat U, V;
     arma::vec s;
     arma::svd(U, s, V, M_c);
     
     // Extract the first left singular vector as the new a
     arma::vec a = U.col(0);
     
     // Update c via soft thresholding using M * a (which should be (p x 1))
     c = ust(M * a, eta);
     
     // Calculate discrepancy and update iteration count
     dis = arma::as_scalar(arma::max(arma::abs(c - c_old)));
     c_old = c;
     ++i;
   }
 } else if(kappa > 0 && kappa < 0.5) {  // Solve for 0 < kappa < 0.5
   double kappa2 = (1 - kappa) / (1 - 2 * kappa);
   
   auto h = [&](double lambda) {
     arma::vec alpha = arma::solve(M + lambda * arma::eye<arma::mat>(p, p), M * c);
     return arma::dot(alpha, alpha) - 1.0 / (kappa2 * kappa2);
   };
   
   while (h(eps) * h(1e30) > 0) {
     while (h(eps) <= 1e5) {
       M *= 2;
       c *= 2;
     }
   }
   
   while (dis > eps && i <= maxstep) {
     if (h(eps) * h(1e30) > 0) {
       while (h(eps) <= 1e5) {
         M *= 2;
         c *= 2;
       }
     }
     
     double lambda = eps;
     for (double l = eps; l < 1e30; l *= 2) {
       if (h(l) * h(eps) < 0) {
         lambda = l;
         break;
       }
     }
     
     arma::vec a = kappa2 * arma::solve(M + lambda * arma::eye<arma::mat>(p, p), M * c);
     c = ust(M * a, eta);
     
     dis = arma::as_scalar(arma::max(arma::abs(c - c_old)));
     c_old = c;
     ++i;
   }
 }
 
 return c;
}