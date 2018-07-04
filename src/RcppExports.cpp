// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// hamming_ngCMatrix_x_only
Rcpp::IntegerVector hamming_ngCMatrix_x_only(Rcpp::S4 obj);
RcppExport SEXP _spamming_hamming_ngCMatrix_x_only(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(hamming_ngCMatrix_x_only(obj));
    return rcpp_result_gen;
END_RCPP
}
// hamming_ngCMatrix_x_and_y
IntegerMatrix hamming_ngCMatrix_x_and_y(Rcpp::S4 objx, Rcpp::S4 objy);
RcppExport SEXP _spamming_hamming_ngCMatrix_x_and_y(SEXP objxSEXP, SEXP objySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type objx(objxSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type objy(objySEXP);
    rcpp_result_gen = Rcpp::wrap(hamming_ngCMatrix_x_and_y(objx, objy));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spamming_hamming_ngCMatrix_x_only", (DL_FUNC) &_spamming_hamming_ngCMatrix_x_only, 1},
    {"_spamming_hamming_ngCMatrix_x_and_y", (DL_FUNC) &_spamming_hamming_ngCMatrix_x_and_y, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spamming(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
