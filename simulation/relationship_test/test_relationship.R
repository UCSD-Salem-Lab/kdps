source("test_kdps.R")

# for(n_relationship in c(100, 500, 1000, 5000, 10000, 50000, 1e5)){
for(n_relationship in c(1e5)){
  # for(seed in c(561892347, 492357816, 634582197)){
  for(seed in c(492357816, 634582197)){
    for(phenotypic_naive in c(TRUE, FALSE)){
      
      cat(paste0("N Relationship: ", n_relationship, "\n"))
      cat(paste0("Seed: ", seed, "\n"))
      cat(paste0("Phenotypic Naive: ", phenotypic_naive, "\n"))
      
      result = test_kdps(
        n_relationship = n_relationship,
        seed = seed,
        fuzziness = 0,
        phenotype_name = "pheno1",
        prioritize_high = FALSE,
        prioritize_low = FALSE,
        phenotype_rank = c("DISEASED", "HEALTHY"),
        phenotypic_naive = phenotypic_naive
      )
      
      filename = paste0("relationship", n_relationship, "_seed", seed, "_", ifelse(phenotypic_naive, "naive", "notnaive"), ".rda")
      # filename = file.path("relationship_test", filename)
      save(result, file = filename)
      
    }
  }
}
