#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Wide Kernel PLS Algorithm Implementation
//'
//' This function implements the wide kernel PLS algorithm as described by
//' Rannar et al. (1994). It's optimized for datasets with more variables than observations.
//'
//' @param X The predictor matrix
//' @param Y The response matrix
//' @param ncomp Number of components to extract
//' @param center Logical, whether to center X and Y (default: true)
//' @param tol Numeric, tolerance used for determining convergence (default: sqrt(machine_eps))
//' @param maxit Integer, maximum number of iterations (default: 100)
//'
//' @return A list containing only:
//' \item{coefficients}{The regression coefficients matrix}
//' \item{projection}{The projection matrix}
//' \item{Xmeans}{The column means of X (if centered)}
//' \item{Ymeans}{The column means of Y (if centered)}
//'
//' @export
// [[Rcpp::export]]
List widekernelpls_fit(arma::mat X, arma::mat Y, int ncomp, 
                      bool center = true,
                      double tol = 0,
                      int maxit = 100) {
  int nobj = X.n_rows;
  int npred = X.n_cols;
  int nresp = Y.n_cols;
  
  // Matrices for scores
  mat TT(nobj, ncomp, fill::zeros);
  mat U(nobj, ncomp, fill::zeros);
   
  // Identity matrix of size nobj
  mat In = eye<mat>(nobj, nobj);
   
  // Means for centering
  rowvec Xmeans(npred, fill::zeros);
  rowvec Ymeans(nresp, fill::zeros);
   
  // Center variables if requested
  if (center) {
    Xmeans = mean(X, 0);
    for (int j = 0; j < npred; j++) {
      X.col(j) -= Xmeans(j);
      }
    
    Ymeans = mean(Y, 0);
    for (int j = 0; j < nresp; j++) {
      Y.col(j) -= Ymeans(j);
      }
    }
  
  // Calculate initial cross-product matrices
  mat XXt = X * X.t();   // tcrossprod in R
  mat YYt = Y * Y.t();   // tcrossprod in R
   
  // Main algorithm loop
  for (int a = 0; a < ncomp; a++) {
    mat XXtYYt = XXt * YYt;
    // Avoid problems with negative eigenvalues due to roundoff errors
    XXtYYt = XXtYYt * XXtYYt;
     
    // Initial values
    vec t_a_old = Y.col(0);
    vec t_a(nobj);
     
    // Iterative calculation of scores vector
    int nit = 0;
    while (true) {
      nit++;
       
      t_a = XXtYYt * t_a_old;
      t_a = t_a / sqrt(as_scalar(t_a.t() * t_a));
       
      double converge_check = sum(abs((t_a - t_a_old) / t_a));
      if (std::isnan(converge_check)) {
        // Handle potential division by zero
        converge_check = sum(abs(t_a - t_a_old));
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
    vec u_a = YYt * t_a;
    vec utmp = u_a / as_scalar(t_a.t() * u_a);
     
    // Calculate wpw
    double wpw = sqrt(as_scalar(utmp.t() * XXt * utmp));
     
    // Save scores
    TT.col(a) = t_a * wpw;
    U.col(a) = utmp * wpw;
     
    // Deflate X and Y cross-product matrices
    mat G = In - t_a * t_a.t();
    XXt = G * XXt * G;
    YYt = G * YYt * G;
    }
  
  // Calculate weights, loadings and projection matrix
  mat W = X.t() * U;
   
  // Normalize weights
  for (int j = 0; j < ncomp; j++) {
    W.col(j) = W.col(j) / sqrt(sum(square(W.col(j))));
    }
  
  // Calculate TT^T(T^TT)^-1 efficiently
  mat TTtTinv(nobj, ncomp);
  for (int j = 0; j < ncomp; j++) {
    TTtTinv.col(j) = TT.col(j) / as_scalar(TT.col(j).t() * TT.col(j));
    }
  
  mat P = X.t() * TTtTinv;
  mat Q = Y.t() * TTtTinv;
   
  // Calculate rotation matrix R
  mat R;
  if (ncomp == 1) {
    // For 1 component, R == W
    R = W;
    } else {
      mat PW = P.t() * W;
      mat PWinv;
      
      if (nresp == 1) {
        // For single-response models, direct calculation of (P^tW)^-1
        PWinv = eye<mat>(ncomp, ncomp);
        vec bidiag(ncomp - 1);
        
        for (int j = 0; j < ncomp - 1; j++) {
          bidiag(j) = -PW(j, j+1);
          }
        
        for (int a = 0; a < ncomp - 1; a++) {
          double cumulative = 1.0;
          for (int j = a; j < ncomp - 1; j++) {
            cumulative *= bidiag(j);
            PWinv(a, j+1) = cumulative;
            }
          }
        } else {
          // Use solve for the general case (equivalent to backsolve in R)
          PWinv = solve(trimatu(PW), eye<mat>(ncomp, ncomp));
          }
        R = W * PWinv;
        }
    
    // Calculate regression coefficients for each number of components
    // B is a 3D array with dimensions [npred, nresp, ncomp]
    cube B(npred, nresp, ncomp);
    
    for (int a = 0; a < ncomp; a++) {
      // Use only the first a+1 components
      B.slice(a) = R.cols(0, a) * Q.cols(0, a).t();
      }
    
    // Return only the requested outputs
    return List::create(
      Named("coefficients") = B,
      Named("projection") = R,
      Named("Xmeans") = Xmeans,
      Named("Ymeans") = Ymeans
    );
    }