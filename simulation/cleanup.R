library(dplyr)
library(writexl)

fl = list.files(path = "fuzziness_test/", pattern = "rda")
all_result_summary = vector()
for(f in fl){
  load(file.path("fuzziness_test/", f))
  d = result
  rm(result)
  
  comp_time = d$runtime
  original_set = d$phenotype_summary$original_set
  unrelated_set = d$phenotype_summary$unrelated_set
  n_remove = d$n_remove
  total_n = d$total_n
  disease_percent = original_set[1]/sum(original_set)
  disease_percent_unrelated = unrelated_set[1]/sum(unrelated_set)
  phenotypic_naive = d$phenotypic_naive
    
  result_summary = tibble(
    fuzziness = gsub(pattern = "fuzziness", replacement = "", unlist(strsplit(f, split = "_"))[1]),
    comp_time = comp_time,
    total_n = total_n,
    n_remove = n_remove,
    disease_percent = disease_percent,
    disease_percent_unrelated = disease_percent_unrelated,
    phenotypic_naive = phenotypic_naive
  )
  
  all_result_summary = rbind.data.frame(all_result_summary, result_summary)
  
}
fuzziness_result_summary = all_result_summary
rm(all_result_summary)

fl = list.files(path = "phenotype_test/", pattern = "rda")
all_result_summary = vector()
for(f in fl){
  load(file.path("phenotype_test/", f))
  d = result
  rm(result)
  
  comp_time = d$runtime
  # original_set = d$phenotype_summary$original_set
  # unrelated_set = d$phenotype_summary$unrelated_set
  n_remove = d$n_remove
  total_n = d$total_n
  # disease_percent = original_set[1]/sum(original_set)
  # disease_percent_unrelated = unrelated_set[1]/sum(unrelated_set)
  phenotypic_naive = d$phenotypic_naive
  
  result_summary = tibble(
    phenotype = gsub(pattern = "phenotype", replacement = "", unlist(strsplit(f, split = "_"))[1]),
    comp_time = comp_time,
    total_n = total_n,
    n_remove = n_remove,
    # disease_percent = disease_percent,
    # disease_percent_unrelated = disease_percent_unrelated,
    phenotypic_naive = phenotypic_naive
  )
  
  all_result_summary = rbind.data.frame(all_result_summary, result_summary)
  
}
phenotype_result_summary = all_result_summary
rm(all_result_summary)

fl = list.files(path = "relationship_test/", pattern = "rda")
all_result_summary = vector()
for(f in fl){
  load(file.path("relationship_test/", f))
  d = result
  rm(result)
  
  comp_time = d$runtime
  original_set = d$phenotype_summary$original_set
  unrelated_set = d$phenotype_summary$unrelated_set
  n_remove = d$n_remove
  total_n = d$total_n
  disease_percent = original_set[1]/sum(original_set)
  disease_percent_unrelated = unrelated_set[1]/sum(unrelated_set)
  phenotypic_naive = d$phenotypic_naive
  
  result_summary = tibble(
    relationship = gsub(pattern = "relationship", replacement = "", unlist(strsplit(f, split = "_"))[1]),
    comp_time = comp_time,
    total_n = total_n,
    n_remove = n_remove,
    disease_percent = disease_percent,
    disease_percent_unrelated = disease_percent_unrelated,
    phenotypic_naive = phenotypic_naive
  )
  
  all_result_summary = rbind.data.frame(all_result_summary, result_summary)
  
}
relationship_result_summary = all_result_summary
rm(all_result_summary)


write_xlsx(fuzziness_result_summary, path = "fuzziness_result_summary.xlsx")
write_xlsx(phenotype_result_summary, path = "phenotype_result_summary.xlsx")
write_xlsx(relationship_result_summary, path = "relationship_result_summary.xlsx")


