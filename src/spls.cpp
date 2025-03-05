#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Declare external functions (assumed to be implemented separately)
Rcpp::List correctp(arma::mat X, arma::mat Y, double eta, int K, double kappa, std::string select, std::string fit);
arma::vec spls_dv(const arma::mat& Z, double eta, double kappa, double eps, int maxstep);
List widekernelpls_fit(arma::mat X, arma::mat Y, int ncomp, bool center = true, double tol = 1.5e-8, int maxit = 100);

// [[Rcpp::export]]
List spls_cpp(arma::mat X, arma::mat Y, int K, double eta, double kappa = 0.5,
              std::string select = "pls2", std::string fit = "simpls",
              bool scale_x = true, bool scale_y = false,
              double eps = 1e-4, int maxstep = 100, bool trace = false) {

  int n = X.n_rows;
  int p = X.n_cols;
  int q = Y.n_cols;
  arma::rowvec meanX, normX, meanY, normY;

  // Center and scale Y
  meanY = arma::mean(Y, 0);
  Y.each_row() -= meanY;
  if (scale_y) {
    normY = arma::stddev(Y, 0, 0);
    Y.each_row() /= normY;
  } else {
    normY.ones(q);
  }

  // Center and scale X
  meanX = arma::mean(X, 0);
  X.each_row() -= meanX;
  if (scale_x) {
    normX = arma::stddev(X, 0, 0);
    if (any(normX < std::numeric_limits<double>::epsilon())) {
      stop("Some of the predictor columns have zero variance.");
    }
    X.each_row() /= normX;
  } else {
    normX.ones(p);
  }

  // Initialize output matrices
  arma::mat betahat(p, q, fill::zeros);
  std::vector<arma::mat> betamat;
  arma::mat X1 = X;
  arma::mat Y1 = Y;

  // Correct input parameters
  List type = correctp(X, Y, eta, K, kappa, select, fit);
  eta = type["eta"];
  K = type["K"];
  kappa = type["kappa"];
  select = Rcpp::as<std::string>(type["select"]);
  fit = Rcpp::as<std::string>(type["fit"]);

  // We also keep track of new predictors in each iteration,
  // but we convert the indices from 0-indexing to 1-indexing.
  std::vector<std::vector<int>> new2As(K);
  arma::uvec ip = arma::linspace<arma::uvec>(0, p - 1, p);

  // Also store the final active set (cumulative A) in an arma::uvec
  arma::uvec finalActiveSet;

  for (int k = 0; k < K; ++k) {

    // Compute Z = X'Y
    arma::mat Z = X1.t() * Y1;

    // Compute sparse direction vector
    arma::vec what = spls_dv(Z, eta, kappa, eps, maxstep);

    // Identify selected predictors; update finalActiveSet each iteration
    finalActiveSet = find((what != 0) || (betahat.col(0) != 0));
    arma::uvec new2A = find((what != 0) && (betahat.col(0) == 0));

    // Convert new2A to a vector of ints and add 1 (to switch from 0-indexing to 1-indexing)
    std::vector<int> tempNew2 = arma::conv_to<std::vector<int>>::from(new2A);
    for (size_t i = 0; i < tempNew2.size(); i++) {
      tempNew2[i] += 1;
    }
    new2As[k] = tempNew2;

    // Fit PLS with selected predictors.
    // Note: finalActiveSet is 0-indexed; when using it to extract rows,
    // it is correct, but when reporting indices to R we add 1.
    arma::mat XA = X.cols(finalActiveSet);
    List plsfit = widekernelpls_fit(XA, Y, std::min(k + 1, (int)finalActiveSet.n_elem), false);

    // Extract coefficients from the PLS fit.
    if (plsfit.containsElementNamed("coefficients")) {
      if (Rf_isMatrix(plsfit["coefficients"])) {
        betahat.zeros();
        betahat.rows(finalActiveSet) = as<arma::mat>(plsfit["coefficients"]);
      } else {
        arma::cube coefCube = as<arma::cube>(plsfit["coefficients"]);
        betahat.zeros();
        betahat.rows(finalActiveSet) = coefCube.slice(0);
      }
    } else {
      stop("Error: 'coefficients' missing in PLS output.");
    }

    betamat.push_back(betahat);

    // Extract projection matrix and update Y1 or X1 based on the selection method.
    if (plsfit.containsElementNamed("projection")) {
      arma::mat pj = as<arma::mat>(plsfit["projection"]);
      if (select == "pls2") {
        Y1 = Y - X * betahat;
      } else if (select == "simpls") {
        arma::mat pw = pj * inv(pj.t() * pj) * pj.t();
        X1.cols(finalActiveSet) -= X.cols(finalActiveSet) * pw;
      }
    } else {
      stop("Error: 'projection' missing in PLS output.");
    }

    // Optional trace output: print new predictors (already converted to 1-indexing)
    if (trace) {
      Rcout << "- " << k + 1 << "th step (K=" << k + 1 << "):\n";
      for (size_t i = 0; i < new2As[k].size() && i < 10; i++) {
        Rcout << new2As[k][i] << " ";
      }
      Rcout << "\n";
    }
  }

  // Return the final active set as a 1-indexed numeric vector (to match the R version)
  arma::vec A_final = conv_to<arma::vec>::from(finalActiveSet) + 1;

  return List::create(
    Named("x") = X,
    Named("y") = Y,
    Named("betahat") = betahat,
    Named("A") = A_final,
    Named("betamat") = betamat,
    Named("new2As") = new2As,
    Named("mu") = meanY,
    Named("meanx") = meanX,
    Named("normx") = normX,
    Named("normy") = normY,
    Named("eta") = eta,
    Named("K") = K,
    Named("kappa") = kappa,
    Named("select") = select,
    Named("fit") = fit
  );
}
