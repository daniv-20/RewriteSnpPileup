// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// run_snp_pileup_logic
void run_snp_pileup_logic(const std::vector<std::string>& input_args);
RcppExport SEXP _snp_plp_run_snp_pileup_logic(SEXP input_argsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type input_args(input_argsSEXP);
    run_snp_pileup_logic(input_args);
    return R_NilValue;
END_RCPP
}
// htslib_version
void htslib_version();
RcppExport SEXP _snp_plp_htslib_version() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    htslib_version();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_snp_plp_run_snp_pileup_logic", (DL_FUNC) &_snp_plp_run_snp_pileup_logic, 1},
    {"_snp_plp_htslib_version", (DL_FUNC) &_snp_plp_htslib_version, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_snp_plp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
