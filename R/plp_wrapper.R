#' SNP Pileup Analysis
#'
#' Runs the SNP pileup analysis with the given arguments.
#'
#' @param file Input file path (e.g., VCF file).
#' @param output Output file path.
#' @param bamfiles Vector of bamfile paths. 
#' @return None. Results are written to the output file.
#' @export
run_snp_pileup <- function(vcffile, ## vcf file path
  output, ## output file path
  bamfiles, ## vector of file paths
  count_orphans = FALSE,
gzipped = FALSE,
ignore_overlaps = FALSE,
min_base_quality = 0,
min_map_quality = 0,
min_read_counts = 0,
max_depth = 4000,
progress = FALSE,
psuedo_snps = 0,
verbose = FALSE,
debug_mode = FALSE
) {
    if (!file.exists(vcffile)) stop("Input file does not exist: ", vcffile)
    ##if (threads <= 0) stop("Threads must be greater than 0.")
  
    qual_args = c("-Q", min_base_quality, 
    "-q", min_map_quality, 
    "-r", min_read_counts, 
    "-d", max_depth,
  "-P", psuedo_snps)
    
  if (count_orphans){
    qual_args = c(qual_args, "-A")
  }
  if(gzipped){
    qual_args = c(qual_args, "-g")
  }
  if(ignore_overlaps){
    qual_args = c(qual_args, "-x")
  }
  if(progress){
    qual_args = c(qual_args, "-p")
  }
  if(verbose){
    qual_args = c(qual_args, "-v")
  }
  if(debug){
    qual_args = c(qual_args, "-d")
  }
    
    args <- c(qual_args, vcffile, output, bamfiles)
    .Call("_snp_plp_run_snp_pileup_logic", args)
}
