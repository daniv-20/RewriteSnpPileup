#devtools::document()

#devtools::install_github("daniv-20/RewriteSnpPileup", ref = "gh_build")

#setwd("/home/nfs/vaithid1/FACETS/RewriteSnpPileup")

BiocManager::install("Rhtslib", force = TRUE)



vcf = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz"




system.file(package = "snp.plp")

# system.file(package = "snp.plp")
# [1] "/home/nfs/vaithid1/R/x86_64-pc-linux-gnu-library/4.3/snp.plp"

rm(list = ls()) 

snp.plp::htsvers()

# snp.plp::test_pl(vcffile = "/scratch/dani/RewriteSnpPileup/data", 
# bamfiles = bams, output = "please_work_finally_thanks", count_orphans = FALSE)

datapath = "/ebfs/epibio/seshanv/snp-pileup"
#vcf = "/scratch/dani/RewriteSnpPileup/data/clean.vcf.gz"
vcf = "/scratch/dani/snpplpextras/data/clean.vcf.gz"

bams = c(file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"))

sink("log.txt")

snp.plp::run_snp_pileup(vcffile = vcf, 
bamfiles = bams, output = "/scratch/dani/snpplpextras/data/bgztest.csv.gz", count_orphans = TRUE, min_read_counts = 5)

remove.packages("snp.plp")

rm(list = ls()) 
.Call("run_snp_pileup_logic")

.Call("htslib_version_cpp")

## find the location of the R header files for the compiler

system("R CMD config --cppflags")
system("R CMD config --ldflags")


system("Rscript -e 'Rhtslib::pkgconfig(\"PKG_LIBS\")'")
## '/opt/R/4.3.1/lib/R/site-library/Rhtslib/usrlib/libhts.a' -lcurl

system.file(package="Rhtslib")
#"/opt/R/4.3.1/lib/R/site-library/Rhtslib"


## Result of above heeader search
# > system("R CMD config --cppflags")
# -I/opt/R/4.3.1/lib/R/include
# > system("R CMD config --ldflags")
# -Wl,--export-dynamic -fopenmp -L/usr/local/lib -L/opt/R/4.3.1/lib/R/lib -lR -lpcre2-8 -llzma -lbz2 -lz -ltirpc -lrt -ldl -lm -licuuc -licui18n