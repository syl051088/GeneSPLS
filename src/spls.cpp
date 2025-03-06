#include "spls.h"
#include "pls.h"
#include "spls_dv.h"
#include "correctp.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name spls_cpp
//' @title Sparse Partial Least Squares Regression
//'
//' @description C++ implementation of the Sparse Partial Least Squares algorithm
//' 
//' @param x Predictor matrix (n x p)
//' @param y Response matrix (n x q)
//' @param K Number of latent components
//' @param eta Sparsity parameter, 0 < eta < 1
//' @param kappa If y is multivariate, 0 < kappa <= 0.5
//' @param select Method for variable selection: "pls2" (update Y) or "simpls" (update X)
//' @param fit Always "widekernelpls" in this implementation
//' @param scale_x Whether to scale the predictor matrix
//' @param scale_y Whether to scale the response matrix
//' @param eps Convergence criterion for direction vector calculation
//' @param maxstep Maximum number of iterations for direction vector calculation
//' @param trace Whether to print progress information
//' 
//' @return A list containing model information including:
//' \item{betahat}{Regression coefficients}
//' \item{A}{Active set of predictors}
//' \item{projection}{Projection matrix}
//' 
//' @useDynLib GeneSPLS, .registration = TRUE
//' @import Rcpp
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List spls_cpp(arma::mat x, arma::mat y, int K, double eta, double kappa = 0.5,
                    std::string select = "pls2", std::string fit = "widekernelpls",
                    bool scale_x = true, bool scale_y = false, 
                    double eps = 1e-4, int maxstep = 100, bool trace = false) {
 
 // Initialize
 int n = x.n_rows;
 int p = x.n_cols;
 int q = y.n_cols;
 arma::rowvec one = arma::ones<arma::rowvec>(n);
 
 // Center & scale y & x
 arma::rowvec mu = one * y / n;
 for (int j = 0; j < q; j++) {
   y.col(j) -= arma::as_scalar(mu(j));
 }
 
 arma::rowvec meanx = one * x / n;
 for (int j = 0; j < p; j++) {
   x.col(j) -= arma::as_scalar(meanx(j));
 }
 
 arma::rowvec normx(p, arma::fill::ones);
 if (scale_x) {
   for (int j = 0; j < p; j++) {
     normx(j) = std::sqrt(arma::sum(arma::square(x.col(j))) / (n-1));
     if (normx(j) < arma::datum::eps) {
       Rcpp::stop("Some of the columns of the predictor matrix have zero variance.");
     }
     x.col(j) /= normx(j);
   }
 }
 
 arma::rowvec normy(q, arma::fill::ones);
 if (scale_y) {
   for (int j = 0; j < q; j++) {
     normy(j) = std::sqrt(arma::sum(arma::square(y.col(j))) / (n-1));
     if (normy(j) < arma::datum::eps) {
       Rcpp::stop("Some of the columns of the response matrix have zero variance.");
     }
     y.col(j) /= normy(j);
   }
 }
 
 // Initialize objects
 arma::mat betahat = arma::zeros<arma::mat>(p, q);
 Rcpp::List betamat(K);
 arma::mat x1 = x;
 arma::mat y1 = y;
 
 // Validate parameters using correctp
 Rcpp::List type = correctp(x, y, eta, K, kappa, select, fit);
 eta = Rcpp::as<double>(type["eta"]);
 K = Rcpp::as<int>(type["K"]);
 kappa = Rcpp::as<double>(type["kappa"]);
 select = Rcpp::as<std::string>(type["select"]);
 fit = "widekernelpls";  // Always use widekernelpls
 
 // Main iteration
 Rcpp::List new2As(K);
 
 if (trace) {
   Rcpp::Rcout << "The variables that join the set of selected variables at each step:\n";
 }
 
 arma::uvec A;
 arma::mat pj;
 
 for (int k = 0; k < K; k++) {
   // Define Z
   arma::mat Z = x1.t() * y1;
   
   // Fit direction vector
   arma::vec what = spls_dv(Z, eta, kappa, eps, maxstep);
   
   // Construct A - indices of non-zero coefficients
   arma::uvec what_nonzero = arma::find(what != 0);
   arma::uvec beta_nonzero = arma::find(betahat.col(0) != 0);
   
   // A is the union of indices where what is non-zero or beta is non-zero
   A = arma::unique(arma::join_cols(what_nonzero, beta_nonzero));
   
   // new2A is indices where what is non-zero and beta is zero
   arma::uvec new2A = arma::find((what != 0) % (betahat.col(0) == 0));
   
   // Fit PLS with predictors in A
   if (A.n_elem > 0) {
     arma::mat xA = x.cols(A);
     int ncomp = std::min(k+1, static_cast<int>(A.n_elem));
     
     // Use widekernelpls_fit (assuming center=false like in the R code)
     Rcpp::List plsfit = widekernelpls_fit(xA, y, ncomp, false, 1.5e-8, 100);
     
     // Extract coefficients for the last component
     arma::cube B = Rcpp::as<arma::cube>(plsfit["coefficients"]);
     arma::mat coeffs = B.slice(ncomp-1);  // Get coefficients for the last component (0-indexed)
     
     // Update betahat
     betahat.zeros();
     for (unsigned int i = 0; i < A.n_elem; i++) {
       betahat.row(A(i)) = coeffs.row(i);
     }
     
     // Get projection matrix
     pj = Rcpp::as<arma::mat>(plsfit["projection"]);
   } else {
     betahat.zeros();
     // Define a default projection matrix if needed
     pj = arma::mat(p, 1, arma::fill::zeros);
   }
   
   // Store in betamat - need to copy betahat to a new matrix for storage
   arma::mat beta_copy = betahat;
   betamat[k] = Rcpp::wrap(beta_copy);
   
   // Update based on selection method
   if (select == "pls2") {
     y1 = y - x * betahat;
   } else if (select == "simpls") {
     if (A.n_elem > 0) {
       arma::mat pw = pj * arma::solve(pj.t() * pj, pj.t());
       x1 = x;
       x1.cols(A) = x.cols(A) - x.cols(A) * pw;
     } else {
       x1 = x;
     }
   }
   
   // Store new variables
   new2As[k] = Rcpp::wrap(new2A + 1);  // Convert to 1-based indexing for R
   
   // Print out variables that join the active set
   if (trace) {
     Rcpp::Rcout << "- " << k+1 << "th step (K=" << k+1 << "):\n";
     for (unsigned int i = 0; i < new2A.n_elem; i++) {
       Rcpp::Rcout << (new2A(i) + 1) << " ";  // Convert to 1-based indexing for display
       if ((i+1) % 10 == 0) Rcpp::Rcout << "\n";
     }
     if (new2A.n_elem > 0) Rcpp::Rcout << "\n";
   }
 }
 
 // Convert A to 1-based indexing for R
 arma::uvec A_1based;
 if (!A.empty()) {
   A_1based = A + 1;
 }
 
 // Return results
 return Rcpp::List::create(
   Rcpp::Named("x") = x,
   Rcpp::Named("y") = y,
   Rcpp::Named("betahat") = betahat,
   Rcpp::Named("A") = Rcpp::wrap(A_1based),
   Rcpp::Named("betamat") = betamat,
   Rcpp::Named("new2As") = new2As,
   Rcpp::Named("mu") = mu,
   Rcpp::Named("meanx") = meanx,
   Rcpp::Named("normx") = normx,
   Rcpp::Named("normy") = normy,
   Rcpp::Named("eta") = eta,
   Rcpp::Named("K") = K,
   Rcpp::Named("kappa") = kappa,
   Rcpp::Named("select") = select,
   Rcpp::Named("fit") = fit,
   Rcpp::Named("projection") = pj
 );
}