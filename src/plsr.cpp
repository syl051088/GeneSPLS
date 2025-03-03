#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Wide Kernel PLS Algorithm Implementation
//'
//' This function implements the wide kernel PLS algorithm as described by
//' Rannar et al. (1994). It's optimized for datasets with more variables than observations.
//'
//' @param X Predictor matrix (n x p), p is the number of predictors
//' @param Y Response matrix (n x m), m is the number of responses
//' @param ncomp Number of components to extract
//' @param center Logical, whether to center X and Y (default: true)
//' @param tol Numeric, tolerance used for determining convergence (default: 1e-6)
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
//' sim_data <- data.frame(Y = Y, X)
//' 
//' pls_cpp = GeneSPLS::widekernelpls_fit(X, Y, 5, center = FALSE)
//' pls_cpp$projection
//' pls_cpp$coefficients
//'
//' @export
// [[Rcpp::export]]
List widekernelpls_fit(arma::mat X,
                       arma::mat Y,
                       int ncomp,
                       bool center = true,
                       double tol = 1.5e-8,
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
    Ymeans = mean(Y, 0);
    for (int j = 0; j < npred; ++j) X.col(j) -= Xmeans(j);
    for (int j = 0; j < nresp; ++j) Y.col(j) -= Ymeans(j);
    }
  
  // Calculate initial cross-product matrices
  mat XXt = X * X.t();
  mat YYt = Y * Y.t();
  
  for (int a = 0; a < ncomp; ++a) {
    
    // Avoid problems with negative eigenvalues due to roundoff errors
    mat XXtYYt = XXt * YYt;
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
          warning("No convergence in %d iterations", maxit);
          break;
          }
        }
    
    // Calculate u vector
    vec u_a = YYt * t_a;
    vec utmp = u_a / as_scalar(t_a.t() * u_a);
    
    double wpw = sqrt(as_scalar(utmp.t() * XXt * utmp));
     
    // Save scores
    TT.col(a) = t_a * wpw;
    U.col(a) = utmp * wpw;
     
    // Deflate X and Y cross-product matrices
    mat G = In - t_a * t_a.t();
    XXt = G * XXt * G;
    YYt = G * YYt * G;
    }
  
  // Calculate weights
  mat W = X.t() * U;
  for (int j = 0; j < ncomp; ++j)
    W.col(j) = W.col(j) / sqrt(sum(square(W.col(j)))); // Normalization
  
  // Calculate TT^T(T^TT)^-1
  mat TTtTinv(nobj, ncomp);
  for (int j = 0; j < ncomp; ++j)
    TTtTinv.col(j) = TT.col(j) / as_scalar(TT.col(j).t() * TT.col(j));
  
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
        
        // For single-response models, directly calculate of (P^tW)^-1
        PWinv = eye<mat>(ncomp, ncomp);
        vec bidiag(ncomp - 1);
        
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
          PWinv = solve(trimatu(PW), eye<mat>(ncomp, ncomp));
          }
        R = W * PWinv;
        }
    
    // Calculate regression coefficients for each number of components
    // B is a 3D array with dimensions [npred, nresp, ncomp]
    cube B(npred, nresp, ncomp);
    
    for (int a = 0; a < ncomp; ++a) {
      // Use only the first a+1 components
      B.slice(a) = R.cols(0, a) * Q.cols(0, a).t();
      }
    
    return List::create(
      Named("coefficients") = B,
      Named("projection") = R,
      Named("Xmeans") = Xmeans,
      Named("Ymeans") = Ymeans
    );
  }

// [OLD] Wide Kernel Partial Least Squares (PLS) Implementation
List widekernelpls_rcpp(arma::mat X, arma::mat Y, int ncomp) {
  
  int n = X.n_rows; // Number of observations
  int p = X.n_cols; // Number of predictors
  int m = Y.n_cols; // Number of response variables
  
  // Ensure ncomp does not exceed available dimensions
  int ncomp_actual = std::min({n, p, ncomp});
  
  // Ensure Y is a column matrix
  if (Y.n_cols == 1) {
    Y.reshape(n, 1);
  }
  
  // Debug: Print dimensions
  Rcpp::Rcout << "X dimensions: " << X.n_rows << "x" << X.n_cols << std::endl;
  Rcpp::Rcout << "Y dimensions: " << Y.n_rows << "x" << Y.n_cols << std::endl;
  Rcpp::Rcout << "Using " << ncomp_actual << " components (min of n, p, and ncomp)" << std::endl;
  
  // Center X and Y
  arma::rowvec X_means = mean(X, 0);
  arma::rowvec Y_means = mean(Y, 0);
  X.each_row() -= X_means;
  Y.each_row() -= Y_means;
  
  // Compute Gram Matrices
  arma::mat XXt = X * X.t();
  arma::mat YYt = Y * Y.t();
  
  // Initialize matrices
  arma::mat T(n, ncomp_actual, fill::zeros); // X scores (Fixed size)
  arma::mat U(n, ncomp_actual, fill::zeros); // Y scores (Fixed size)
  arma::mat B(p, m, fill::zeros);     // Regression coefficients
  arma::mat W(p, ncomp_actual, fill::zeros); // Loading weights
  arma::mat P(p, ncomp_actual, fill::zeros); // Loadings
  arma::mat Q(ncomp_actual, m, fill::zeros); // Y loadings
  
  // Debugging print statements
  Rcpp::Rcout << "T dimensions: " << T.n_rows << "x" << T.n_cols << std::endl;
  Rcpp::Rcout << "U dimensions: " << U.n_rows << "x" << U.n_cols << std::endl;
  Rcpp::Rcout << "W dimensions: " << W.n_rows << "x" << W.n_cols << std::endl;
  Rcpp::Rcout << "P dimensions: " << P.n_rows << "x" << P.n_cols << std::endl;
  
  // PLS Component Iterations
  for (int a = 0; a < ncomp_actual; ++a) {
    Rcpp::Rcout << "Processing component " << a << " / " << ncomp_actual - 1 << std::endl;
    Rcpp::Rcout << "T size: " << T.n_rows << "x" << T.n_cols << ", Accessing T.col(" << a << ")" << std::endl;
    Rcpp::Rcout << "U size: " << U.n_rows << "x" << U.n_cols << ", Accessing U.col(" << a << ")" << std::endl;
    
    arma::mat XXtYYt = XXt * YYt;
    XXtYYt = XXtYYt * XXtYYt;
    
    arma::vec t_a = Y.col(0);
    t_a.reshape(n, 1);  // Ensure correct shape
    
    arma::vec t_a_old = t_a;
    double tol = 1e-6;
    int max_iter = 100, iter = 0;
    
    // Power method to compute leading eigenvector
    while (iter < max_iter) {
      iter++;
      t_a = XXtYYt * t_a_old;
      t_a /= norm(t_a, 2);
      if (norm(t_a - t_a_old, 2) < tol) break;
      t_a_old = t_a;
    }
    
    arma::vec u_a = YYt * t_a;
    u_a.reshape(n, 1);  // Ensure correct shape
    
    double norm_factor = dot(t_a, u_a);
    arma::vec u_tmp = u_a / norm_factor;
    
    double wpw = sqrt(dot(u_tmp, XXt * u_tmp));
    T.col(a) = t_a * wpw;
    U.col(a) = u_tmp * wpw;
    
    // Deflate XXt and YYt
    arma::mat G = eye(n, n) - t_a * t_a.t();
    XXt = G * XXt * G;
    YYt = G * YYt * G;
  }
  
  // Return results
  return List::create(
    Named("coefficients") = B,
    Named("scores") = T,
    Named("loadings") = P,
    Named("projection") = W
  );
}
