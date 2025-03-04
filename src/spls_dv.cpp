#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Univariate Soft Thresholding (UST) Function
arma::vec ust(const arma::vec& b, double eta) {
  arma::vec b_ust = arma::zeros<arma::vec>(b.n_elem);
  if (eta < 1) {
    arma::vec valb = abs(b) - eta * max(abs(b));
    uvec idx = find(valb >= 0);
    b_ust.elem(idx) = valb.elem(idx) % sign(b.elem(idx));
  }
  return b_ust;
}

// Sparse PLS Direction Vector Estimation
// [[Rcpp::export]]
arma::vec spls_dv(const arma::mat& Z, double eta, double kappa, double eps, int maxstep) {
  int p = Z.n_rows;
  int q = Z.n_cols;
  double Znorm1 = median(abs(vectorise(Z)));
  arma::mat Z_scaled = Z / Znorm1;

  // Case: Univariate response
  if(q == 1) {
    return ust(Z_scaled.col(0), eta);
  }

  // Case: Multivariate response
  arma::mat M = Z_scaled * Z_scaled.t();  // (p x p)

  // Initialize c as in the original R code: (p x 1) vector with constant value 10
  arma::vec c = 10 * ones<vec>(p);
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
      svd(U, s, V, M_c);

      // Extract the first left singular vector as the new a
      arma::vec a = U.col(0);

      // Update c via soft thresholding using M * a (which should be (p x 1))
      c = ust(M * a, eta);

      // Calculate discrepancy and update iteration count
      dis = as_scalar(max(abs(c - c_old)));
      c_old = c;
      i++;
    }
  } else if(kappa > 0 && kappa < 0.5) {  // Solve for 0 < kappa < 0.5
    double kappa2 = (1 - kappa) / (1 - 2 * kappa);

    auto h = [&](double lambda) {
      arma::vec alpha = solve(M + lambda * eye<mat>(p, p), M * c);
      return dot(alpha, alpha) - 1.0 / (kappa2 * kappa2);
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

      arma::vec a = kappa2 * solve(M + lambda * eye<mat>(p, p), M * c);
      c = ust(M * a, eta);

      dis = as_scalar(max(abs(c - c_old)));
      c_old = c;
      i++;
    }
  }

  return c;
}
