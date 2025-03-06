#include "pls.h"
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name widekernelpls_fit
//' @title Wide Kernel PLS Algorithm Implementation
//'
//' @description This function implements the wide kernel PLS algorithm as described by
//' Rannar et al. (1994). It's optimized for datasets with more variables than observations.
//'
//' @param X Predictor matrix (n x p), p is the number of predictors
//' @param Y Response matrix (n x m), m is the number of responses
//' @param ncomp Number of components to extract
//' @param center Logical, whether to center X and Y (default: true)
//' @param tol Numeric, tolerance used for determining convergence (default: 1.5e-8)
//' @param maxit Integer, maximum number of iterations (default: 100)
//'
//' @return A list containing:
//' \item{coefficients}{The regression coefficients matrix}
//' \item{projection}{The projection matrix}
//' \item{Xmeans}{The column means of X (if centered)}
//' \item{Ymeans}{The column means of Y (if centered)}
//' 
//' @useDynLib GeneSPLS, .registration = TRUE
//' @import Rcpp
//' 
//' @examples
//' 
//' set.seed(815)
//' 
//' # Simulating high-dimensional predictor matrix (p >> n)
//' n <- 10   # Number of observations
//' p <- 20  # Number of predictors
//' 
//' X <- matrix(rnorm(n * p), nrow = n, ncol = p)  # Random normal predictors
//' 
//' # True coefficients: Only first 5 predictors contribute to Y
//' beta <- c(rnorm(5, sd = 2), rep(0, p - 5))
//' 
//' # Generating response variable Y with N(0,1) noise
//' Y <- X %*% beta + rnorm(n)
//' 
//' 
//' pls_cpp = GeneSPLS::widekernelpls_fit(X, Y, 5, center = FALSE)
//' pls_cpp$projection
//' pls_cpp$coefficients
//'
//' @export
// [[Rcpp::export]]
Rcpp::List widekernelpls_fit(arma::mat X,
                             arma::mat Y,
                             int ncomp,
                             bool center = true,
                             double tol = 1.5e-8,
                             int maxit = 100) {
 int nobj = X.n_rows;
 int npred = X.n_cols;
 int nresp = Y.n_cols;
 
 // Matrices for scores
 arma::mat TT(nobj, ncomp, arma::fill::zeros);
 arma::mat U(nobj, ncomp, arma::fill::zeros);
 
 // Identity matrix of size nobj
 arma::mat In = arma::eye<arma::mat>(nobj, nobj);
 
 // Means for centering
 arma::rowvec Xmeans(npred, arma::fill::zeros);
 arma::rowvec Ymeans(nresp, arma::fill::zeros);
 
 // Center variables if requested
 if (center) {
   Xmeans = arma::mean(X, 0);
   Ymeans = arma::mean(Y, 0);
   for (int j = 0; j < npred; ++j) X.col(j) -= Xmeans(j);
   for (int j = 0; j < nresp; ++j) Y.col(j) -= Ymeans(j);
 }
 
 // Calculate initial cross-product matrices
 arma::mat XXt = X * X.t();
 arma::mat YYt = Y * Y.t();
 
 for (int a = 0; a < ncomp; ++a) {
   
   // Avoid problems with negative eigenvalues due to roundoff errors
   arma::mat XXtYYt = XXt * YYt;
   XXtYYt = XXtYYt * XXtYYt;
   
   // Initial values
   arma::vec t_a_old = Y.col(0);
   arma::vec t_a(nobj);
   
   // Iterative calculation of scores vector
   int nit = 0;
   while (true) {
     nit++;
     
     t_a = XXtYYt * t_a_old;
     t_a = t_a / std::sqrt(arma::as_scalar(t_a.t() * t_a));
     
     double converge_check = arma::sum(arma::abs((t_a - t_a_old) / t_a));
     if (std::isnan(converge_check)) {
       // Handle potential division by zero
       converge_check = arma::sum(arma::abs(t_a - t_a_old));
     }
     
     if (converge_check < tol) {
       break;
     } else {
       t_a_old = t_a;
     }
     
     if (nit >= maxit) {
       Rcpp::warning("No convergence in %d iterations", maxit);
       break;
     }
   }
   
   // Calculate u vector
   arma::vec u_a = YYt * t_a;
   arma::vec utmp = u_a / arma::as_scalar(t_a.t() * u_a);
   
   double wpw = std::sqrt(arma::as_scalar(utmp.t() * XXt * utmp));
   
   // Save scores
   TT.col(a) = t_a * wpw;
   U.col(a) = utmp * wpw;
   
   // Deflate X and Y cross-product matrices
   arma::mat G = In - t_a * t_a.t();
   XXt = G * XXt * G;
   YYt = G * YYt * G;
 }
 
 // Calculate weights
 arma::mat W = X.t() * U;
 for (int j = 0; j < ncomp; ++j)
   W.col(j) = W.col(j) / std::sqrt(arma::sum(arma::square(W.col(j)))); // Normalization
 
 // Calculate TT^T(T^TT)^-1
 arma::mat TTtTinv(nobj, ncomp);
 for (int j = 0; j < ncomp; ++j)
   TTtTinv.col(j) = TT.col(j) / arma::as_scalar(TT.col(j).t() * TT.col(j));
 
 arma::mat P = X.t() * TTtTinv;
 arma::mat Q = Y.t() * TTtTinv;
 
 // Calculate rotation matrix R
 arma::mat R;
 if (ncomp == 1) {
   // For 1 component, R == W
   R = W;
 } else {
   arma::mat PW = P.t() * W;
   arma::mat PWinv;
   
   if (nresp == 1) {
     // For single-response models, directly calculate of (P^tW)^-1
     PWinv = arma::eye<arma::mat>(ncomp, ncomp);
     arma::vec bidiag(ncomp - 1);
     
     for (int j = 0; j < ncomp - 1; ++j) bidiag(j) = -PW(j, j+1);
     
     for (int a = 0; a < ncomp - 1; ++a) {
       double cumulative = 1.0;
       
       for (int j = a; j < ncomp - 1; ++j) {
         cumulative *= bidiag(j);
         PWinv(a, j+1) = cumulative;
       }
     }
   } else {
     // Use solve for the general case
     PWinv = arma::solve(arma::trimatu(PW), arma::eye<arma::mat>(ncomp, ncomp));
   }
   R = W * PWinv;
 }
 
 // Calculate regression coefficients for each number of components
 // B is a 3D array with dimensions [npred, nresp, ncomp]
 arma::cube B(npred, nresp, ncomp);
 
 for (int a = 0; a < ncomp; ++a) {
   // Use only the first a+1 components
   B.slice(a) = R.cols(0, a) * Q.cols(0, a).t();
 }
 
 return Rcpp::List::create(
   Rcpp::Named("coefficients") = B,
   Rcpp::Named("projection") = R,
   Rcpp::Named("Xmeans") = Xmeans,
   Rcpp::Named("Ymeans") = Ymeans
 );
}