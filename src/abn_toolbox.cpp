#include <Rcpp.h>
#include <math.h>

//' @title logit and logit functions
//' 
//' @description transform \code{x} either via the logit, or expit.
//'
//'
//' @param x a numeric vector
//' @export
//' @rdname logit
// [[Rcpp::export]]
Rcpp::NumericVector logit_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector result(n);

  for(int i = 0; i < n; ++i) { 
    result[i] = log( x[i] / (1.0 - x[i]) );
  }

  return result;
}

//' @export
//' @rdname logit
// [[Rcpp::export]]
Rcpp::NumericVector expit_cpp(Rcpp::NumericVector x) { 
  int n = x.size();
  Rcpp::NumericVector result(n);

  for (int i=0; i < n; ++i) { 
    result[i] = 1.0 / (1.0 + exp (-1.0 * x[i]));
  }
  return result;
}
