// #include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include "snp-pileup-rev.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib> // For exit()
#include <cstdio>  // For printf, FILE
#include <sstream>
#include <ctime>

// using namespace Rcpp;

// Define the arguments structure
#ifndef SNP_PILEUP_REV_H_DEFINED
#define SNP_PILEUP_REV_H_DEFINED
struct arguments
{
  std::vector<std::string> args;
  bool count_orphans = true;
  bool ignore_overlaps = false;
  int min_base_quality = 0;
  int min_map_quality = 0;
  std::vector<int> min_read_counts;
  int max_depth = 4000;
  int rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
  BGZF *gzippedPointer = nullptr;
  bool debug_mode = false;
};
#endif

void gzip_output_dbg(arguments arguments, std::string str, FILE *fp) {
    if (!arguments.gzippedPointer) {
        std::cerr << "Error: BGZF pointer is NULL!" << std::endl;
        return;
    }
    
    int written = bgzf_write(arguments.gzippedPointer, str.c_str(), str.length());
    if (written < 0) {
        std::cerr << "Error: bgzf_write failed!" << std::endl;
    } else {
        std::cout << "Successfully wrote " << written << " bytes to gzip file" << std::endl;
    }
    bgzf_flush(arguments.gzippedPointer);
}

// Checks if a string ends with a certain substring
inline bool ends_with(std::string const &value, std::string const &ending)
{
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

//' @export
//' @export
extern "C" int dbg_gzip(arguments arguments)
{
    std::cout << "In Main yay" << std::endl;
    clock_t start = clock();

    // Check if there are enough arguments
    if (arguments.args.size() < 2) {
        std::cerr << "Error: Expected at least 2 arguments, but received " << arguments.args.size() << std::endl;
        return 1;
    }

    std::string fname = arguments.args[1];  // Now safe to access

    if (fname.empty()) {
        std::cerr << "Error: Output filename is empty!" << std::endl;
        return 1;
    }

    if (!ends_with(fname, ".gz")) {
        fname += ".gz";
    }

    // Open gzip file
    arguments.gzippedPointer = bgzf_open(fname.c_str(), "w+");
    if (!arguments.gzippedPointer) {
        std::cerr << "Error: Failed to open gzipped output file: " << fname << std::endl;
        return 1;
    }
    std::cout << "Opened gzipped file: " << fname << std::endl;

    std::ostringstream output;
    output << "Chromosome,Position,Ref,Alt\n";
    std::string header = output.str();

    // Use valid FILE* instead of nullptr
    gzip_output_dbg(arguments, header, arguments.gzippedPointer);

    // Simulate writing output data
    for (int i = 0; i < 5; i++) {
        std::ostringstream row;
        row << "chr1," << i+1 << ",A,T\n";
        gzip_output_dbg(arguments, row.str(), arguments.gzippedPointer);
    }

    // Flush and close the file safely
    int flush_status = bgzf_flush(arguments.gzippedPointer);
    if (flush_status < 0) {
        std::cerr << "Error: bgzf_flush() failed!" << std::endl;
    }

    bgzf_close(arguments.gzippedPointer);
    return 0;
}
