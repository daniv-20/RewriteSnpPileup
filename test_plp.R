#devtools::document()

#devtools::install_github("daniv-20/RewriteSnpPileup", ref = "gh_build")

#setwd("/home/nfs/vaithid1/FACETS/RewriteSnpPileup")

datapath = "/ebfs/epibio/seshanv/snp-pileup"

vcf = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz"
bams = c(file.path(datapath, "HCC1143_BC10.bam"), file.path(datapath, "HCC1143_BL10.bam"))




rm(list = ls()) 

snp.plp::htsvers()

snp.plp::test_pl

snp.plp::run_snp_pileup(vcffile = "/usr/local/share/VCF/common_all_20180418_dedup_bg.vcf.gz", 
bamfiles = c("bam", "bam bam"), output = "please_work_finally_thanks", count_orphans = 1.5)



remove.packages("snp.plp")

rm(list = ls()) 
.Call("run_snp_pileup_logic")

.Call("htslib_version_cpp")