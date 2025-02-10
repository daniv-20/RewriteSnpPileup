#' SNP Pileup Analysis
#'
#' Runs the SNP pileup analysis with the given arguments.
#'
#' @param vcffile Input file path (e.g., VCF file).
#' @param output Output file path.
#' @param bamfiles Vector of bamfile paths.
#' @param  count_orphans Do not discard anomalous read pairs.
#' @param max_depth Sets the maximum depth. Default is 4000.
##gzipped Compresses the output file with BGZF. If false, output file is CSV.
#' @param psuedo_snps Every n positions, if there is no SNP, insert a blank record with the total count at the position.
#' @param min_map_quality Sets the minimum threshold for mapping quality. The default is 0.
#' @param min_base_quality Sets the minimum threshold for base quality. Default is 0.
#' @param min_read_counts Comma separated list of minimum read counts for position to be output. Default is 0.
#' @param ignore_overlaps If true, disables read-pair overlap detection.
#' @param debug_mode If true, shows many detailed messages.
#' @return None. Results are written to the output file.
#' @export
run_snp_pileup <- function(
  vcffile, ## vcf file path
  output, ## output file path
  bamfiles, ## vector of file paths
  count_orphans = TRUE,
  # gzipped = FALSE,
  ignore_overlaps = FALSE,
  min_base_quality = 0,
  min_map_quality = 0,
  min_read_counts = 0,
  max_depth = 4000,
  # progress = FALSE,
  psuedo_snps = 0,
  # verbose = FALSE,
  debug_mode = FALSE) {

  print("Entering snp.plp function")

  if(length(min_read_counts) < length(bamfiles)){
    if(length(min_read_counts) == 1){
      count = min_read_counts
      min_read_counts = rep(count, length(bamfiles))
    }
    else{
      min_read_counts = c(min_read_counts, rep(0, length(bamfiles) - length(min_read_counts)))
    }
  }
  if(length(min_read_counts) > length(bamfiles)){
    min_read_counts = min_read_counts[1:length(bamfiles)]
  }
  
args <- list(
  count_orphans = ifelse(count_orphans, 1L, 0L),
  ignore_overlaps =  ifelse(ignore_overlaps, 1L, 0L),
  min_base_quality = as.integer(min_base_quality),
  min_map_quality = as.integer(min_map_quality),
  min_read_counts = as.integer(min_read_counts),
  max_depth = as.integer(max_depth),
  psuedo_snps = as.integer(psuedo_snps),
  args = c(vcffile, output, bamfiles)
)

print("Made args list")
print(paste0(names(args), collapse = ", "))

.Call("run_snp_pileup_logic", args) ## will change to _facets_run_snp_pileup_logic
  
}

#' @export
htsvers = function(){
  print("Calling htslib version")
  .Call("htslib_version_cpp")
}