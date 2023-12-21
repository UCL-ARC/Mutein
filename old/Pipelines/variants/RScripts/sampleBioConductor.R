if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

# install necessary libraries
if (!require("GenomicRanges", quietly = TRUE))
    BiocManager::install("GenomicRanges")
if (!require("VariantAnnotation", quietly = TRUE))
    BiocManager::install("VariantAnnotation")
if (!require("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
    BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
if (!require("PolyPhen.Hsapiens.dbSNP131", quietly = TRUE))
    BiocManager::install("PolyPhen.Hsapiens.dbSNP131")


# https://bioconductor.org/packages/release/workflows/vignettes/sequencing/inst/doc/sequencing.html

library(GenomicRanges)
GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 3, 4)),
      IRanges(1:10, width=5), strand='-',
      score=101:110, GC = runif(10))


#https://bioconductor.org/packages/release/workflows/vignettes/variants/inst/doc/Annotating_Genomic_Variants.html#amino-acid-coding-changes-in-non-synonymous-variants


# SETUP
library(VariantAnnotation)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(PolyPhen.Hsapiens.dbSNP131)

# Exploring variants
#hardcoded_vcf_path = "/home/rachel/UCL/github/MuteinData/vcf_keogh/keogh2018.vcf"
print(getwd())
setwd("\\\\wsl$\\Ubuntu\\home\\rachel\\UCL\\github\\MuteinData\\vcf_keogh")
setwd("C:/UCL/temp_r")
print(getwd())
filenames <- list.files(pattern="*.*", full.names=TRUE)
print(filenames)

hardcoded_vcf_path = "C:/UCL/temp_r/keogh2018.vcf"
#file <- system.file("vcf", hardcoded_vcf_path)
vcf <-system.file("vcf",hardcoded_vcf_path)
#vcf<-readVcf(hardcoded_vcf_path,"hg19")
hdr <- scanVcfHeader(vcf)
print(info(hdr))