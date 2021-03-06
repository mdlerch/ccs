// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// larscpp
arma::mat larscpp(arma::mat x, arma::mat beta, arma::vec y, int maxk);
RcppExport SEXP ccs_larscpp(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP maxkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type maxk(maxkSEXP);
    rcpp_result_gen = Rcpp::wrap(larscpp(x, beta, y, maxk));
    return rcpp_result_gen;
END_RCPP
}
// lassocpp
arma::mat lassocpp(arma::mat x, arma::mat beta, arma::vec y, int maxk);
RcppExport SEXP ccs_lassocpp(SEXP xSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP maxkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type maxk(maxkSEXP);
    rcpp_result_gen = Rcpp::wrap(lassocpp(x, beta, y, maxk));
    return rcpp_result_gen;
END_RCPP
}
