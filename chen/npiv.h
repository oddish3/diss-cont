// chen/npiv.h
#ifndef NPIV_H
#define NPIV_H

#include <RcppArmadillo.h>

// Declare the npiv function
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y);

#endif // NPIV_H
