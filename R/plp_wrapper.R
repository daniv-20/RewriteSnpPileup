#' SNP Pileup Analysis
#'
#' Runs the SNP pileup analysis with the given arguments.
#'
#' @param vcffile Input file path (e.g., VCF file).
#' @param output Output file path.
#' @param bamfiles Vector of bamfile paths.
#' @param count_orphans If true, do not discard anomalous read pairs. Default TRUE.
#' @param ignore_overlaps If true, disables read-pair overlap detection. Default FALSE.
#' @param max_depth Sets the maximum depth. Default is 4000.
#' @param psuedo_snps Every n positions, if there is no SNP, insert a blank record with the total count at the position. Default is 0.
#' @param min_map_quality Sets the minimum threshold for mapping quality. The default is 0.
#' @param min_base_quality Sets the minimum threshold for base quality. Default is 0.
#' @param min_read_counts Comma separated list of minimum read counts for position to be output. Default is 0.
#' @return None. Results are written to the output file.
#' @export
run_snp_pileup <- function(
  vcffile, ## vcf file path
  output, ## output file path
  bamfiles, ## vector of file paths
  count_orphans = TRUE,
  ignore_overlaps = FALSE,
  max_depth = 4000,
  psuedo_snps = 0,
  min_map_quality = 0,
  min_base_quality = 0,
  min_read_counts = 0
  ) {

  ## make sure we have the correct number of min_read_counts

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
.Call("run_snp_pileup_logic", args)   
}

#' @title htsvers function
#' @description Gets the version of HTSLIB being run
#' @return Returns the version of htslib being run. 
#' @export
htsvers = function(){
  # print("Calling htslib version")
  .Call("htslib_version_cpp")
}