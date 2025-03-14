\name{run_snp_pileup}
\alias{run_snp_pileup}
\title{SNP Pileup Analysis}
\description{
Runs the SNP pileup analysis with the given arguments.
}
\usage{
run_snp_pileup(
  vcffile,
  output,
  bamfiles,
  count_orphans = TRUE,
  ignore_overlaps = FALSE,
  max_depth = 4000,
  psuedo_snps = 0,
  min_map_quality = 0,
  min_base_quality = 0,
  min_read_counts = 0
)
}
\arguments{
  \item{vcffile}{Input file path (e.g., VCF file).}
  \item{output}{Output file path.}
  \item{bamfiles}{Vector of BAM file paths.}
  \item{count_orphans}{If true, do not discard anomalous read pairs. Default \code{TRUE}.}
  \item{ignore_overlaps}{If true, disables read-pair overlap detection. Default \code{FALSE}.}
  \item{max_depth}{Sets the maximum depth. Default is \code{4000}.}
  \item{psuedo_snps}{Every n positions, if there is no SNP, insert a blank record with the total count at the position. Default is \code{0}.}
  \item{min_map_quality}{Sets the minimum threshold for mapping quality. Default is \code{0}.}
  \item{min_base_quality}{Sets the minimum threshold for base quality. Default is \code{0}.}
  \item{min_read_counts}{Comma-separated list of minimum read counts for position to be output. Default is \code{0}.}
}
\value{
None. Results are written to the output file.
}
\examples{
run_snp_pileup(
  vcffile = "example.vcf",
  output = "output.txt",
  bamfiles = c("sample1.bam", "sample2.bam"),
  count_orphans = TRUE,
  ignore_overlaps = FALSE,
  max_depth = 4000,
  psuedo_snps = 0,
  min_map_quality = 20,
  min_base_quality = 30,
  min_read_counts = c(5, 10)
)
}
\keyword{SNP}
