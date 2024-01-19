source("test_kdps.R")

for(fuzziness in c(0,1,2,5,10)){
  for(seed in c(561892347, 492357816, 634582197)){
    for(phenotypic_naive in c(TRUE, FALSE)){
      
      cat(paste0("Fuzziness: ", fuzziness, "\n"))
      cat(paste0("Seed: ", seed, "\n"))
      cat(paste0("Phenotypic Naive: ", phenotypic_naive, "\n"))
      
      result = test_kdps(
        n_relationship = 10000,
        seed = seed,
        fuzziness = fuzziness,
        phenotype_name = "pheno1",
        prioritize_high = FALSE,
        prioritize_low = FALSE,
        phenotype_rank = c("DISEASED", "HEALTHY"),
        phenotypic_naive = phenotypic_naive
      )
      
      filename = paste0("fuzziness", fuzziness, "_seed", seed, "_", ifelse(phenotypic_naive, "naive", "notnaive"), ".rda")
      filename = file.path("fuzziness_test", filename)
      save(result, file = filename)
      
    }
  }
}
