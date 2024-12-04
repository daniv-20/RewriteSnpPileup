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

outfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "outputs", "test1.csv")

input_args = c(vcf, outfile, file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"))


# Compile the file using sourceCpp
Rcpp::sourceCpp("bgzip_snp-pileup.cpp", rebuild = TRUE, verbose = TRUE)

htslib_version()

delete_if_exists(outfile)

run_snp_pileup(input_args)
# Cleanup: Restore the original Makevars
Sys.unsetenv("R_MAKEVARS_USER")
unlink(temp_makevars)

