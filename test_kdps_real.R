library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(progress)
library(knitr)

source("kdps.R")

phenotype_file = "data/real.sample"
kinship_file = "data/real.kin0"

# d = as.data.frame(fread(phenotype_file))
# fid = strsplit(d$FID_IID, split = "[_]") %>%
#   lapply(function(x){x[1]}) %>%
#   unlist()
# iid = strsplit(d$FID_IID, split = "[_]") %>%
#   lapply(function(x){x[2]}) %>%
#   unlist()
# 
# d = data.frame(
#   FID = fid,
#   IID = iid,
#   dbgapid = d$dbgapid,
#   DIAB = d$DIAB
# )
# 
# fwrite(d, file = "real.sample", sep = "\t")

test_real = kdps(phenotype_file = "data/real.sample", # Real data .gitignored
                 kinship_file = "data/real.kin0", # Real data .gitignored
                 fuzziness = 0,
                 phenotype_name = "DIAB",
                 prioritize_high = FALSE,
                 prioritize_low = FALSE,
                 phenotype_rank = c(1, 0),
                 # Header names in the phenotype file
                 fid_name = "FID",
                 iid_name = "IID",
                 # Header names in the kinship file
                 fid1_name = "#FID1",
                 iid1_name = "IID1",
                 fid2_name = "FID2",
                 iid2_name = "IID2",
                 kinship_name = "KINSHIP",
                 kinship_threshold = 0.0442,
                 phenotypic_naive = FALSE
)
