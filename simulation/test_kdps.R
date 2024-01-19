library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(progress)
library(knitr)

kinship = fread("../data/kinship.txt")
pheno = fread("../data/pheno.txt")
source("../kdps.R")

test_kdps = function(
    n_relationship = 10000,
    seed = 492357816,
    fuzziness = 0,
    phenotype_name = "pheno2",
    prioritize_high = FALSE,
    prioritize_low = FALSE,
    phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY"),
    phenotypic_naive = FALSE
){
  
  set.seed(seed)
  
  kinship_subset = kinship[sample(1:dim(kinship)[1], n_relationship), ]
  pheno_subest = pheno %>%
    filter(IID %in% unique(c(kinship_subset$IID1, kinship_subset$IID2)))
  
  kinship_table = kinship_subset %>%
    mutate(fid1_iid1 = paste(FID1, IID1, sep = "_")) %>%
    mutate(fid2_iid2 = paste(FID2, IID2, sep = "_")) %>%
    mutate(fid1_iid1 = paste0("subject_", fid1_iid1)) %>%
    mutate(fid2_iid2 = paste0("subject_", fid2_iid2)) %>%
    mutate(related = KINSHIP >= 0.0442) %>%
    filter(related) %>%
    select(fid1_iid1, fid2_iid2, related)
  
  relationship_count = table(c(
    kinship_table[["fid1_iid1"]],
    kinship_table[["fid2_iid2"]]
  )) %>%
    table()
  
  relationship_count = tibble(
    n_relatedness = names(relationship_count),
    count = as.vector(relationship_count)
  )
  
  print(kable(relationship_count))
  
  fwrite(kinship_subset, file = "kinship_temp.txt", sep = " ")
  fwrite(pheno_subest, file = "pheno_temp.txt", sep = " ")
  
  start_time = Sys.time()
  subject_to_remove = kdps(phenotype_file = "pheno_temp.txt", 
                           kinship_file = "kinship_temp.txt",
                           fuzziness = fuzziness,
                           phenotype_name = phenotype_name,
                           prioritize_high = prioritize_high,
                           prioritize_low = prioritize_low,
                           phenotype_rank = phenotype_rank, 
                           phenotypic_naive = phenotypic_naive)
  end_time = Sys.time()
  runtime = end_time - start_time
  cat(paste0("Total runtime: ", round(runtime, 3), " seconds.\n"))
  
  unrelated_phenotype = pheno_subest %>%
    filter(!(IID %in% subject_to_remove[["IID"]])) %>%
    select(FID, IID, all_of(phenotype_name))
  
  if(is.numeric(unlist(unrelated_phenotype[,3]))){
    phenotype_summary = list(
      original_set = summary(pheno_subest[[phenotype_name]]),
      unrelated_set = summary(unrelated_phenotype[,3])
    )
  }else{
    phenotype_summary = list(
      original_set = table(pheno_subest[[phenotype_name]]),
      unrelated_set = table(unrelated_phenotype[,3])
    )
  }
  
  test_result = list(
    n_relationship = n_relationship,
    seed = seed,
    fuzziness = fuzziness,
    phenotype_name = phenotype_name,
    prioritize_high = prioritize_high,
    prioritize_low = prioritize_low,
    phenotype_rank = phenotype_rank,
    phenotypic_naive = phenotypic_naive,
    runtime = runtime,
    total_n = dim(pheno_subest)[1],
    n_remove = dim(subject_to_remove)[1],
    phenotype_summary = phenotype_summary
  )
  
  file.remove("kinship_temp.txt")
  file.remove("pheno_temp.txt")
  
  return(test_result)
}
