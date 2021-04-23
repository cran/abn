#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;

// [[Rcpp::export]]
double factorial(double n)
{
  return (n == 1.0 || n == 0.0) ? 1 : factorial(n - 1.0) * n;
}


// [[Rcpp::export]]

Rcpp::List irls_poisson_cpp(arma::mat A, arma::vec b, double maxit, double tol)
{
//Def
arma::vec x;
x.zeros(A.n_cols,1);
arma::vec xold;
arma::mat varmatrix;

double nobs;
nobs = A.n_rows;
double ll;
arma::vec e;
double ssr;
double aic;
double bic;
double df;
double mdl;


arma::vec W(nobs);
arma::vec unit(nobs);
unit.ones(nobs);
arma::vec eta(nobs);
arma::vec g(nobs);
arma::vec f(nobs);
arma::vec gprime(nobs);
arma::vec z(nobs);

for (int i = 0; i < maxit; ++i) {
  eta = A * x;
  
//  g = exp(eta);
//  gprime = exp(eta);
   for (int j=0; j < nobs; ++j) {
     g[j] = exp(eta[j]);
//     gprime[j] = exp(eta[j]);
   }
    gprime = g;
    
  z = eta+(b-g)/gprime;
  
  W = gprime % gprime;
  W /= g;
  xold = x;
  
  //coefficients
  varmatrix = A.t()*(W % A.each_col());
  x = arma::solve(varmatrix, A.t()*(W % z), arma::solve_opts::no_approx);
  
if(sqrt(arma::dot(x-xold,x-xold)) < tol){
 break;
}}

df = A.n_cols;

//loglik
    
    for (int j = 0; j < nobs; ++j) {
      f[j] = log(factorial(1.0 * b[j]));
      }
    ll = arma::accu(b % (A * x) - exp(A * x) - f);
    
aic = - 2 * ll + 2 * df;

bic = - 2 * ll + log(nobs) * df;
    
//sse
e = (b - A*x);
    //TEST GK 21.03.2021
//ssr = accu(e.t()*e);
    ssr = arma::accu(e.t()*e);
    
    
mdl = 1;
    
//return
return Rcpp::List::create(
  Rcpp::Named("coefficients") = x,
  Rcpp::Named("loglik") = ll,
  Rcpp::Named("aic") = aic,
  Rcpp::Named("bic") = bic,
  Rcpp::Named("mdl") = mdl,
  Rcpp::Named("sse") = ssr,
  Rcpp::Named("varcov") = varmatrix
);
}



