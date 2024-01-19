library(dplyr)
library(data.table)

set.seed(492357816)

kinship = fread("ukb_rel.dat")
names(kinship) = c("IID1", "IID2", "HetHet", "IBS0", "KINSHIP")

kinship = kinship %>%
  mutate(FID1 = 0) %>%
  mutate(FID2 = 0) %>%
  select(FID1, IID1, FID2, IID2, HetHet, IBS0, KINSHIP) %>%
  filter(IID1 > 0) %>%
  filter(IID2 > 0)

samples = c(paste(kinship$FID1, kinship$IID1, sep = "_"),
            paste(kinship$FID2, kinship$IID2, sep = "_")) %>%
  unique()
n_sample = length(samples)

pheno = tibble(
  FID = 0,
  IID = strsplit(samples, split = "[_]") %>%
    lapply(function(x){x[2]}) %>%
    unlist()
)

pheno = pheno %>%
  mutate(pheno1 = sample(c(
    rep("DISEASED",round(0.2*n_sample,0)),
    rep("HEALTHY",n_sample-round(0.2*n_sample,0))
  ))) %>%
  mutate(pheno2 = sample(c(
    rep("DISEASED1",round(0.1*n_sample,0)),
    rep("DISEASED2",round(0.2*n_sample,0)),
    rep("HEALTHY",n_sample-round(0.1*n_sample,0)-round(0.2*n_sample,0)
    )))) %>%
  mutate(pheno3 = round(rnorm(n_sample, mean = 100, sd = 10),2))


fwrite(kinship, file = "kinship.txt", sep = " ")
fwrite(pheno, file = "pheno.txt", sep = " ")
