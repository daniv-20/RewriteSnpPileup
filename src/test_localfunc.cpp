#include <R.h>
#include <Rinternals.h>
#include <iostream>
#include <vector>
#include "zlib.h"
#include "zconf.h"
#include <Rdefines.h>
// #include <Rcpp.h>
// using namespace Rcpp;

struct arguments // some of these things will need to leave here...
{
    bool count_orphans;
    int min_base_quality;
    int min_map_quality;
    std::vector<int> min_read_counts;
    int max_depth;
    int pseudo_snps;
    //BGZF *gzippedPointer;  // Pointer to a BGZF file (part of HTSlib)
    bool gzipped;          // Boolean to indicate gzipped status
    std::vector<std::string> args;  // Store additional arguments as strings
};


extern "C" SEXP localfunc() {
    Rprintf("Hello from localfunc\n");
    return R_NilValue;
}

// extern "C" SEXP extvers() {
//     return Rf_mkString(hts_version());
// }

extern "C" SEXP zlibVers() {
	return Rf_mkString(zlibVersion());
}


// extern "C" SEXP pass_through(SEXP arguments_sexp){

//     List arguments(arguments_sexp);

//     CharacterVector args = arguments["args"];
//     Rcpp::Rcout << "Received named list with elements:\n";
    
//     // Optionally: Return the list back unchanged
//     return arguments_sexp;

// }
//' @export
extern "C" SEXP process_list(SEXP input_list) {
    Rprintf("Hello from process_list\n");

    if (!Rf_isNewList(input_list)) {
        Rf_error("Input must be a named list.");
    }

    // Populate the struct from the SEXP list
    arguments args;
    args.count_orphans = INTEGER(VECTOR_ELT(input_list, 0))[0];
    std::cout << "count_orphans: " << args.count_orphans << "\n";
    args.min_base_quality = INTEGER(VECTOR_ELT(input_list, 1))[0];
    std::cout << "min_base_quality: " << args.min_base_quality << "\n";
    args.min_map_quality = INTEGER(VECTOR_ELT(input_list, 2))[0];
    std::cout << "min_map_quality: " << args.min_map_quality << "\n";

    // Convert SEXP to std::vector<int> for min_read_counts
    SEXP read_counts = VECTOR_ELT(input_list, 3);
    int len = LENGTH(read_counts);
    for (int i = 0; i < len; ++i) {
        args.min_read_counts.push_back(INTEGER(read_counts)[i]);
        std::cout << "min_read_counts[" << i << "]: " << args.min_read_counts[i] << "\n";
    }



    args.max_depth = INTEGER(VECTOR_ELT(input_list, 4))[0];
   std::cout << "max_depth: " << args.max_depth << "\n";

    args.pseudo_snps = INTEGER(VECTOR_ELT(input_list, 5))[0];
    std::cout << "pseudo_snps: " << args.pseudo_snps << "\n";
    // Extract gzipped value
    args.gzipped = LOGICAL(VECTOR_ELT(input_list, 6))[0];
    std::cout << "gzipped: " << args.gzipped << "\n";

    // Convert SEXP to std::vector<std::string> for args
    SEXP char_args = VECTOR_ELT(input_list, 7);
    int char_len = LENGTH(char_args);
    for (int i = 0; i < char_len; ++i) {
        args.args.push_back(CHAR(STRING_ELT(char_args, i)));
        std::cout << "args[" << i << "]: " << args.args[i] << "\n";
    }

    

    // Debug output
    
    
    
   

    // Return the original list for now
    return input_list;
}

// extern "C" SEXP process_list(SEXP input_list){
//     Rprintf("Hello from process_list\n");

//     if (!Rf_isNewList(input_list)) {
//         Rprintf("Input must be a named list.");
//     }

//      // Access elements directly by name
//        // Access elements by index
//     // Populate the struct from the SEXP list
//     arguments args;

//     // Populate the struct from the SEXP list
//     args.count_orphans = LOGICAL(VECTOR_ELT(input_list, 0))[0];       // Convert SEXP to bool
//     args.min_base_quality = INTEGER(VECTOR_ELT(input_list, 1))[0];    // Convert SEXP to int
//     args.min_map_quality = INTEGER(VECTOR_ELT(input_list, 2))[0];     // Convert SEXP to int

//     // Convert SEXP to std::vector<int> for min_read_counts
//     SEXP read_counts = VECTOR_ELT(input_list, 3);
//     int len = LENGTH(read_counts);
//     for (int i = 0; i < len; ++i) {
//         args.min_read_counts.push_back(INTEGER(read_counts)[i]);
//     }

//     args.max_depth = INTEGER(VECTOR_ELT(input_list, 4))[0];           // Convert SEXP to int
//     args.pseudo_snps = INTEGER(VECTOR_ELT(input_list, 5))[0];         // Convert SEXP to int
//     //args.gzippedPointer = nullptr;  // Set pointer to null for now
//     args.gzipped = LOGICAL(VECTOR_ELT(input_list, 6))[0];


//     // Debug output
//     std::cout << "count_orphans: " << args.count_orphans << "\n";
//     std::cout << "gzipped: " << args.gzipped << "\n";
//     std::cout << "min_base_quality: " << args.min_base_quality << "\n";

//     std::cout << "args: ";
//     for (const auto& arg : args.args) {
//         std::cout << arg << " ";
//     }
//     std::cout << "\n";

//     // Return the original list for now
//     return input_list;
// }