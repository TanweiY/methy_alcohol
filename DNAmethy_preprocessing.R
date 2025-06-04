## we used the preprocessed data from our previous project, calculating methylation-based smoking signatures
## The preprocessing methods (including normalization) were the same for both blood and tumor samples
source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

### tumor methylation
tumorm_raw <- loadData(idatPath = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/Validation_cpgs/Data_methylation450")
methy_normILM <- normalizeData(RGset=tumorm_raw, normMethod = "ILM")
save(methy_normILM,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/methy_normILM.RData')

### blood methylation preprocessing
bloodm_raw <- loadData(idatPath = "/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/blood_methylation/rawdata")
bloodm_ILM <- normalizeData(RGset=bloodm_raw, normMethod = "ILM")
save(bloodm_ILM,
     file = '/omics/odcf/analysis/OE0167_projects/dachs_genetic_data_platform/smoking_methylation/processed_data/bloodm_ILM.RData')



