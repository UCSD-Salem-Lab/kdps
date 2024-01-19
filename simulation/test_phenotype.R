source("test_kdps.R")

for(seed in c(561892347, 492357816, 634582197)){
  for(phenotypic_naive in c(TRUE, FALSE)){
    
    cat(paste0("Phenotype: Pheno1", "\n"))
    cat(paste0("Seed: ", seed, "\n"))
    cat(paste0("Phenotypic Naive: ", phenotypic_naive, "\n"))
    
    result = test_kdps(
      n_relationship = 10000,
      seed = seed,
      fuzziness = 0,
      phenotype_name = "pheno1",
      prioritize_high = FALSE,
      prioritize_low = FALSE,
      phenotype_rank = c("DISEASED", "HEALTHY"),
      phenotypic_naive = phenotypic_naive
    )
    
    filename = paste0("phenotypePheno1", "_seed", seed, "_", ifelse(phenotypic_naive, "naive", "notnaive"), ".rda")
    filename = file.path("phenotype_test", filename)
    save(result, file = filename)
    
  }
}

for(seed in c(561892347, 492357816, 634582197)){
  for(phenotypic_naive in c(TRUE, FALSE)){
    
    cat(paste0("Phenotype: Pheno2", "\n"))
    cat(paste0("Seed: ", seed, "\n"))
    cat(paste0("Phenotypic Naive: ", phenotypic_naive, "\n"))
    
    result = test_kdps(
      n_relationship = 10000,
      seed = seed,
      fuzziness = 0,
      phenotype_name = "pheno2",
      prioritize_high = FALSE,
      prioritize_low = FALSE,
      phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY"),
      phenotypic_naive = phenotypic_naive
    )
    
    filename = paste0("phenotypePheno2", "_seed", seed, "_", ifelse(phenotypic_naive, "naive", "notnaive"), ".rda")
    filename = file.path("phenotype_test", filename)
    save(result, file = filename)
    
  }
}

for(seed in c(561892347, 492357816, 634582197)){
  for(phenotypic_naive in c(TRUE, FALSE)){
    
    cat(paste0("Phenotype: Pheno3", "\n"))
    cat(paste0("Seed: ", seed, "\n"))
    cat(paste0("Phenotypic Naive: ", phenotypic_naive, "\n"))
    
    result = test_kdps(
      n_relationship = 10000,
      seed = seed,
      fuzziness = 0,
      phenotype_name = "pheno3",
      prioritize_high = TRUE,
      prioritize_low = FALSE,
      phenotype_rank = c("DISEASED", "HEALTHY"),
      phenotypic_naive = phenotypic_naive
    )
    
    filename = paste0("phenotypePheno3", "_seed", seed, "_", ifelse(phenotypic_naive, "naive", "notnaive"), ".rda")
    filename = file.path("phenotype_test", filename)
    save(result, file = filename)
    
  }
}

