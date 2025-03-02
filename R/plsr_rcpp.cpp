#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Wide Kernel Partial Least Squares (PLS) Implementation
// [[Rcpp::export]]
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
