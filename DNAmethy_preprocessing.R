## we used the preprocessed data from our previous project calculating methylation-based smoking signatures
## the preprocessing methods (including normalization) was the same for both blood and tumor samples
source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")

library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(htmlTable)
library(rmarkdown)
library(devtools)
library(EpiSmokEr)
library(readxl)
library(cgwtools)
library(writexl)

library(minfi)
library(limma)
library(impute)
library(preprocessCore)
library(BiocParallel)

### blood methylation preprocessing
bloodm_raw <- loadData(idatPath = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/blood_methylation/rawdata")
bloodm_ILM <- normalizeData(RGset=bloodm_raw, normMethod = "ILM")
save(bloodm_ILM,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/bloodm_ILM.RData')

### tumor methylation
tumorm_raw <- loadData(idatPath = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/Data_methylation450")
tumorm_ILM <- normalizeData(RGset=tumorm_raw, normMethod = "ILM")
save(tumorm_ILM,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/tumorm_ILM.RData')



