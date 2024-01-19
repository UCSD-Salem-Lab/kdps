library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(progress)
library(knitr)

kinship = fread("simulation/simple_kinship.txt")
pheno = fread("simulation/simple_pheno.txt")
source("kdps.R")

rm(subject_to_remove)
subject_to_remove = kdps(phenotype_file = "simulation/simple_pheno.txt",
                         kinship_file = "simulation/simple_kinship.txt",
                         fuzziness = 0,
                         phenotype_name = "pheno3",
                         prioritize_high = FALSE,
                         prioritize_low = TRUE,
                         phenotype_rank = c("DISEASED", "HEALTHY"),
                         # Header names in the phenotype file
                         fid_name = "FID",
                         iid_name = "IID",
                         # Header names in the kinship file
                         fid1_name = "FID1",
                         iid1_name = "IID1",
                         fid2_name = "FID2",
                         iid2_name = "IID2",
                         kinship_name = "KINSHIP",
                         kinship_threshold = 0.0442,
                         phenotypic_naive = FALSE)
kable(subject_to_remove)
