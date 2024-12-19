devtools::document()

devtools::install_github("daniv-20/RewriteSnpPileup", ref = "gh_build")

setwd("/home/nfs/vaithid1/FACETS/RewriteSnpPileup")

datapath = "/ebfs/epibio/seshanv/snp-pileup"

vcf = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz"
bams = c(file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"))



snp.plp::run_snp_pileup(vcffile = vcf, output = "test_test5", bamfiles = bams, verbose = TRUE, debug_mode = TRUE, progress = TRUE, gzipped = TRUE)
 

