#' @export
test_pl <- function(
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
  
    args <- list(
      args = c(vcffile, output, bamfiles),
      count_orphans = count_orphans, ##ifelse(count_orphans, 1L, 0L),
      ignore_overlaps =  ifelse(ignore_overlaps, 1L, 0L),
      min_base_quality = min_base_quality,
      min_map_quality = min_map_quality,
      min_read_counts = min_read_counts,
      max_depth = max_depth,
      psuedo_snps = psuedo_snps,
      debug_mode =  ifelse(debug_mode, 1L, 0L)
    )
    
    print("Made args list")

    .Call("process_list2", args) ## will change to _facets_run_snp_pileup_logic
  
  }