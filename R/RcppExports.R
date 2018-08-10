# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title logit and logit functions
#' 
#' @description transform \code{x} either via the logit, or expit.
#'
#'
#' @param x a numeric vector
#' @export
#' @rdname logit
logit_cpp <- function(x) {
    .Call('_abn_logit_cpp', PACKAGE = 'abn', x)
}

#' @export
#' @rdname logit
expit_cpp <- function(x) {
    .Call('_abn_expit_cpp', PACKAGE = 'abn', x)
}

irls_binomial_cpp_br <- function(A, b, maxit, tol) {
    .Call('_abn_irls_binomial_cpp_br', PACKAGE = 'abn', A, b, maxit, tol)
}

irls_binomial_cpp_fast_br <- function(A, b, maxit, tol) {
    .Call('_abn_irls_binomial_cpp_fast_br', PACKAGE = 'abn', A, b, maxit, tol)
}

irls_binomial_cpp_fast <- function(A, b, maxit, tol) {
    .Call('_abn_irls_binomial_cpp_fast', PACKAGE = 'abn', A, b, maxit, tol)
}

irls_binomial_cpp <- function(A, b, maxit, tol) {
    .Call('_abn_irls_binomial_cpp', PACKAGE = 'abn', A, b, maxit, tol)
}

irls_gaussian_cpp_fast <- function(A, b, maxit, tol) {
    .Call('_abn_irls_gaussian_cpp_fast', PACKAGE = 'abn', A, b, maxit, tol)
}

irls_gaussian_cpp <- function(A, b, maxit, tol) {
    .Call('_abn_irls_gaussian_cpp', PACKAGE = 'abn', A, b, maxit, tol)
}

factorial_fast <- function(n) {
    .Call('_abn_factorial_fast', PACKAGE = 'abn', n)
}

irls_poisson_cpp_fast <- function(A, b, maxit, tol) {
    .Call('_abn_irls_poisson_cpp_fast', PACKAGE = 'abn', A, b, maxit, tol)
}

factorial <- function(n) {
    .Call('_abn_factorial', PACKAGE = 'abn', n)
}

irls_poisson_cpp <- function(A, b, maxit, tol) {
    .Call('_abn_irls_poisson_cpp', PACKAGE = 'abn', A, b, maxit, tol)
}

mi_cpp <- function(joint_dist) {
    .Call('_abn_mi_cpp', PACKAGE = 'abn', joint_dist)
}
