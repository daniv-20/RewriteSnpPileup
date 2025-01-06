#include <Rcpp.h>
using namespace Rcpp;

// Accept SEXP but use Rcpp::List internally
extern "C" SEXP processNamedListSEXP(SEXP input_sexp) {
    try {
        // Convert SEXP to Rcpp::List
        List input(input_sexp);

        // Access elements by name
        int a = input["a"];
        double b = input["b"];
        std::string c = input["c"];

        // Convert back to SEXP and return
        return wrap(input_sexp);
    } catch (std::exception &ex) {
        Rf_error("Error: %s", ex.what());
    } catch (...) {
        Rf_error("Unknown error");
    }
    return R_NilValue;  // Fallback in case of error
}
