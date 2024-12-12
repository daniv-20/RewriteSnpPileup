rm(list = ls())

getwd()

setwd("/home/nfs/vaithid1/FACETS/RewriteSnpPileup")

source("/home/nfs/vaithid1/FACETS/RewriteSnpPileup/idl.R")

print(system.file(package = "Rhtslib"))


library(snp.plp)

##BiocManager::install("Rhtslib")
#devtools::install_github("daniv-20/RewriteSnpPileup", ref="gh_build")

## get ready for testiiiiing
datapath = "/ebfs/epibio/seshanv/snp-pileup"

vcf = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz"


## set up log file

log = generate_output_file(base_name = "Log", ".txt")

logfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "Logs", log)

## out = generate_output_file("test", "")

outfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "outputs", "ArgsTest_gh_snpplp.csv")

qual_args = c("-d", "2500", "-q", "15", "-Q", "20", "-r", "10")

input_args = c(qual_args, vcf, outfile, file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"), "-p")

##sink(logfile)

cat("Current run: Try progress bar...")
cat(paste("\nArgs: ", paste0(input_args, collapse = ", ")))

# Compile the file using sourceCpp
Rhtslib:::htsVersion()

snp.plp::htslib_version()

delete_if_exists(outfile)

run_snp_pileup(input_args)
# Cleanup: Restore the original Makevars
Sys.unsetenv("R_MAKEVARS_USER")
unlink(temp_makevars)

sink()
cat("unsunk")

## check out results for funsies -> gzipped file

Rhtslib:::htsVersion()
