#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List widekernelpls_fit(arma::mat X,
                       arma::mat Y,
                       int ncomp,
                       bool center = true,
                       std::string validation = "none",
                       int segments = 10, // Default: 10-fold CV
                       std::string segment_type = "random", // "random", "consecutive", "interleaved"
                       double tol = 1.5e-8,
                       int maxit = 100) {

  int nobj = X.n_rows;
  int npred = X.n_cols;
  int nresp = Y.n_cols;

  // Mean centering
  rowvec Xmeans(npred, fill::zeros);
  rowvec Ymeans(nresp, fill::zeros);

  if (center) {
    Xmeans = mean(X, 0);
    Ymeans = mean(Y, 0);
    X.each_row() -= Xmeans;
    Y.each_row() -= Ymeans;
  }

  // If no cross-validation, proceed with the standard wide kernel PLS
  if (validation == "none") {

    mat TT(nobj, ncomp, fill::zeros);
    mat U(nobj, ncomp, fill::zeros);

    mat In = eye<mat>(nobj, nobj);
    mat XXt = X * X.t();
    mat YYt = Y * Y.t();

    for (int a = 0; a < ncomp; ++a) {

      mat XXtYYt = XXt * YYt;
      XXtYYt = XXtYYt * XXtYYt;

      vec t_a_old = Y.col(0);
      vec t_a(nobj);

      int nit = 0;
      while (true) {
        nit++;
        t_a = XXtYYt * t_a_old;
        t_a = t_a / sqrt(as_scalar(t_a.t() * t_a));

        double converge_check = sum(abs((t_a - t_a_old) / t_a));
        if (std::isnan(converge_check)) {
          converge_check = sum(abs(t_a - t_a_old));
        }

        if (converge_check < tol) break;
        else t_a_old = t_a;

        if (nit >= maxit) {
          warning("No convergence in %d iterations", maxit);
          break;
        }
      }

      vec u_a = YYt * t_a;
      vec utmp = u_a / as_scalar(t_a.t() * u_a);
      double wpw = sqrt(as_scalar(utmp.t() * XXt * utmp));

      TT.col(a) = t_a * wpw;
      U.col(a) = utmp * wpw;

      mat G = In - t_a * t_a.t();
      XXt = G * XXt * G;
      YYt = G * YYt * G;
    }

    mat W = X.t() * U;
    for (int j = 0; j < ncomp; ++j)
      W.col(j) = W.col(j) / sqrt(sum(square(W.col(j))));

    mat TTtTinv(nobj, ncomp);
    for (int j = 0; j < ncomp; ++j)
      TTtTinv.col(j) = TT.col(j) / as_scalar(TT.col(j).t() * TT.col(j));

    mat P = X.t() * TTtTinv;
    mat Q = Y.t() * TTtTinv;

    mat R;
    if (ncomp == 1) {
      R = W;
    } else {
      mat PW = P.t() * W;
      mat PWinv = solve(trimatu(PW), eye<mat>(ncomp, ncomp));
      R = W * PWinv;
    }

    cube B(npred, nresp, ncomp);
    for (int a = 0; a < ncomp; ++a) {
      B.slice(a) = R.cols(0, a) * Q.cols(0, a).t();
    }

    return List::create(
      Named("coefficients") = B,
      Named("projection") = R,
      Named("Xmeans") = Xmeans,
      Named("Ymeans") = Ymeans
    );
  }

  // -----------------------------------------------------------
  // CROSS-VALIDATION (mvrCv)
  // -----------------------------------------------------------
  int folds = (validation == "LOO") ? nobj : segments;
  mat PRESS(nresp, ncomp, fill::zeros);
  mat RMSEP(nresp, ncomp, fill::zeros);

  // Generate cross-validation splits
  std::vector<uvec> cv_segments(folds);
  for (int i = 0; i < folds; i++) {
    cv_segments[i] = regspace<uvec>(i, folds, nobj - 1);
  }

  for (int k = 0; k < folds; k++) {

    uvec test_idx = cv_segments[k];
    uvec train_idx = linspace<uvec>(0, nobj - 1, nobj);
    train_idx.shed_rows(test_idx);

    mat X_train = X.rows(train_idx);
    mat Y_train = Y.rows(train_idx);
    mat X_test = X.rows(test_idx);
    mat Y_test = Y.rows(test_idx);

    List model = widekernelpls_fit(X_train, Y_train, ncomp, false, "none", 0, "none");
    cube B_cv = model["coefficients"];

    mat Y_pred(X_test.n_rows, nresp, fill::zeros);
    for (int a = 0; a < ncomp; a++) {
      Y_pred = X_test * B_cv.slice(a);
      mat residuals = Y_test - Y_pred;
      PRESS.col(a) += sum(square(residuals), 0).t();
    }
  }

  RMSEP = sqrt(PRESS / nobj);

  return List::create(
    Named("PRESS") = PRESS,
    Named("RMSEP") = RMSEP
  );
}
