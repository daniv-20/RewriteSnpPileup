rm(list = ls())

getwd()

setwd("/home/nfs/vaithid1/FACETS/RewriteSnpPileup")

source("/home/nfs/vaithid1/FACETS/RewriteSnpPileup/idl.R")

print(system.file(package = "Rhtslib"))

##BiocManager::install("Rhtslib")

temp_makevars <- tempfile()
file.create(temp_makevars)

Sys.setenv(R_MAKEVARS_USER = temp_makevars)

# Set the PKG_CPPFLAGS and PKG_LIBS explicitly for this compilation -> link to rhtslib
Sys.setenv(
 PKG_CPPFLAGS = "-I\"/home/nfs/vaithid1/R/x86_64-pc-linux-gnu-library/4.4/Rhtslib/include\"",
 PKG_LIBS = "-L\"/home/nfs/vaithid1/R/x86_64-pc-linux-gnu-library/4.4/Rhtslib/usrlib\" -lhts"
)

## get ready for testiiiiing
datapath = "/ebfs/epibio/seshanv/snp-pileup"

vcf = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz"


## set up log file

log = generate_output_file(base_name = "Log", ".txt")

logfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "Logs", log)

## out = generate_output_file("test", "")

outfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "outputs", "TestGzipped.csv")

input_args = c(vcf, outfile, file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"), "--debug", "-g")

sink(logfile)

cat("Current run: Try Gzip output, made code update to hopefully fix writer: ")
cat(paste("\nArgs: ", paste0(input_args, collapse = ", ")))

# Compile the file using sourceCpp
Rcpp::sourceCpp("bgzip_snp-pileup.cpp", rebuild = TRUE)

htslib_version()

# delete_if_exists(outfile)

run_snp_pileup(input_args)
# Cleanup: Restore the original Makevars
Sys.unsetenv("R_MAKEVARS_USER")
unlink(temp_makevars)

sink()
cat("unsunk")

## check out results for funsies
##source("/home/nfs/vaithid1/FACETS/RewriteSnpPileup/lookAtOutput.R")
# output = read.csv(outfile)

# ## output columns
# ## chromosome position ref alt  then for each file there is an R E A D

# output.long = pivot_longer(output, cols = !c("Chromosome", "Position", "Ref", "Alt"), names_to = c("file", ".value"),  # `file` becomes one column, RAED stay as values
# names_pattern = "File(\\d+)([RAED])"  # Regex to separate file number and RAED
# ) %>%
#   mutate(across(
#     .cols = -c(Chromosome, Position, Ref, Alt),  # Exclude Ref and Alt columns
#     .fns = as.numeric      # Apply as.numeric
#   )) %>% tbl_summary(
#     by = file,
#     digits = all_continuous() ~ 1,
#     statistic = all_continuous() ~ "{mean}, ({sd})"
#   )

