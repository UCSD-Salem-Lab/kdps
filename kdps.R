library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(progress)

kdps = function(phenotype_file = "data/pheno.txt",
                kinship_file = "data/kinship.txt",
                fuzziness = 0,
                phenotype_name = "pheno2",
                prioritize_high = FALSE,
                prioritize_low = FALSE,
                phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY"),
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
                phenotypic_naive = FALSE){
  
  ### Make phenotype data frame
  pheno = as.data.frame(fread(phenotype_file))
  
  if(!any(prioritize_high, prioritize_low)){
    phenotype = pheno[[phenotype_name]]
    mapping = setNames(length(phenotype_rank):1, phenotype_rank)
    wt = mapping[phenotype]
    names(wt) = NULL
    wt[is.na(wt)] = 0
  }else{
    if(prioritize_high){
      phenotype = as.numeric(pheno[[phenotype_name]])
      na_wt = min(phenotype, na.rm = TRUE) - 1
      wt = ifelse(is.na(phenotype), na_wt, phenotype)
      wt = order(wt)
    }
    if(prioritize_low){
      phenotype = as.numeric(pheno[[phenotype_name]])
      na_wt = max(phenotype, na.rm = TRUE) + 1
      wt = ifelse(is.na(phenotype), na_wt, phenotype)
      wt = order(wt, decreasing = TRUE)
    }
  }
  
  fid_iid = paste0("subject_", paste(pheno[[fid_name]], 
                                     pheno[[iid_name]], 
                                     sep = "_"))
  
  pheno = tibble(
    fid_iid = fid_iid,
    wt = wt
  )
  
  ### Make kinship data frame
  kinship = as.data.frame(fread(kinship_file))
  kinship = tibble(
    fid1 = kinship[[fid1_name]],
    iid1 = kinship[[iid1_name]],
    fid2 = kinship[[fid2_name]],
    iid2 = kinship[[iid2_name]],
    kinship = kinship[[kinship_name]]
  ) 
  
  kinship = kinship %>%
    mutate(fid1_iid1 = paste(fid1, iid1, sep = "_")) %>%
    mutate(fid2_iid2 = paste(fid2, iid2, sep = "_")) %>%
    mutate(fid1_iid1 = paste0("subject_", fid1_iid1)) %>%
    mutate(fid2_iid2 = paste0("subject_", fid2_iid2)) %>%
    mutate(related = kinship >= kinship_threshold) %>%
    filter(related) %>%
    select(fid1_iid1, fid2_iid2, related) 
  
  # Filter out singular relatedness
  cat("Filtering out singular relatedness...\n")
  relationship = table(c(
    kinship[["fid1_iid1"]],
    kinship[["fid2_iid2"]]
  ))
  one_timers = names(relationship[relationship == 1])
  
  singular_nodes = kinship %>%
    filter(fid1_iid1 %in% one_timers) %>%
    filter(fid2_iid2 %in% one_timers) %>%
    left_join(mutate(pheno, fid1_iid1 = fid_iid, wt1 = wt) %>% 
                select(fid1_iid1, wt1),
              by = "fid1_iid1") %>%
    left_join(mutate(pheno, fid2_iid2 = fid_iid, wt2 = wt) %>% 
                select(fid2_iid2, wt2),
              by = "fid2_iid2") %>%
    mutate(subject_to_remove = ifelse(wt2 > wt1, fid1_iid1, fid2_iid2))
  
  if(phenotypic_naive){
    singular_nodes = singular_nodes %>%
      mutate(subject_to_remove = ifelse(runif(dim(singular_nodes)[1]) > 0.5,
                                        fid1_iid1, fid2_iid2))
  }
  
  subject_to_remove_list = singular_nodes[["subject_to_remove"]]
  
  kinship = kinship %>%
    filter(!(fid1_iid1 %in% subject_to_remove_list)) %>%
    filter(!(fid2_iid2 %in% subject_to_remove_list))
  
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  if(dim(kinship)[1] == 0){
    cat("Only singular relatedness are present in the dataset...\n")
    average_wt = (pheno %>%
                    filter(!(fid_iid %in% subject_to_remove_list)))[["wt"]] %>%
      mean()
    cat(paste0("Average subject phenotypic weight: ", round(average_wt,2), "\n\n"))
    
    subject_to_remove_list = tibble(
      FID = strsplit(subject_to_remove_list, split = "_") %>%
        lapply(function(x){x[2]}) %>% unlist(),
      IID = strsplit(subject_to_remove_list, split = "_") %>%
        lapply(function(x){x[3]}) %>% unlist()
    )
    
    return(subject_to_remove_list)
  }
  
  # Filter out singular nodes if fuzziness is 0
  if(fuzziness == 0){
    relationship = table(c(
      kinship[["fid1_iid1"]],
      kinship[["fid2_iid2"]]
    ))
    isolated_n_timer = vector()
    for(n_relation in 2:10){
      cat(paste0("Filtering out isolated super-subjects with ", n_relation, " relatives...\n"))
      n_timers = names(relationship[relationship == n_relation])
      if(length(n_timers) > 0){
        pb = progress_bar$new(total = length(n_timers))
        for(subject in n_timers){
          subject_node = kinship %>%
            filter(fid1_iid1 == subject | fid2_iid2 == subject)
          connected_subjects = c(subject_node[["fid1_iid1"]],
                                 subject_node[["fid2_iid2"]])
          connected_subjects = connected_subjects[connected_subjects != subject]
          if(all(connected_subjects %in% one_timers)){
            isolated_n_timer = rbind.data.frame(
              isolated_n_timer,
              tibble(subject = subject, n = n_relation)
            )
          }
          pb$tick()
        }
      }
    }
    
    subject_to_remove_list = c(subject_to_remove_list, 
                               isolated_n_timer[["subject"]])
    kinship = kinship %>%
      filter(!(fid1_iid1 %in% subject_to_remove_list)) %>%
      filter(!(fid2_iid2 %in% subject_to_remove_list))
  }
  
  # Pruning complex relateness networks
  cat("Pruning complex relateness networks...\n")
  relationship = table(c(
    kinship[["fid1_iid1"]],
    kinship[["fid2_iid2"]]
  ))
  max_count = max(relationship)
  
  while(max_count > 1){
    relationship = table(c(
      kinship[["fid1_iid1"]],
      kinship[["fid2_iid2"]]
    ))
    
    max_count = max(relationship)
    max_count_corrected = max_count - fuzziness
    
    subject_to_remove = (tibble(
      fid_iid = names(relationship[relationship >= max_count_corrected])
    ) %>%
      left_join(pheno, by = "fid_iid") %>%
      arrange(wt))[["fid_iid"]][1]
    
    if(phenotypic_naive){
      subject_to_remove = (tibble(
        fid_iid = names(relationship[relationship >= max_count_corrected])
      ) %>%
        left_join(pheno, by = "fid_iid"))[["fid_iid"]]
      subject_to_remove = subject_to_remove[sample(1:length(subject_to_remove), 1)]
    }
    
    kinship = kinship %>%
      filter(fid1_iid1 != subject_to_remove) %>%
      filter(fid2_iid2 != subject_to_remove)
    
    cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
    subject_to_remove_list = c(subject_to_remove_list, subject_to_remove)
  }
  
  # Work up
  cat("Performing work-up...\n")
  singular_nodes = kinship %>%
    left_join(mutate(pheno, fid1_iid1 = fid_iid, wt1 = wt) %>% 
                select(fid1_iid1, wt1),
              by = "fid1_iid1") %>%
    left_join(mutate(pheno, fid2_iid2 = fid_iid, wt2 = wt) %>% 
                select(fid2_iid2, wt2),
              by = "fid2_iid2") %>%
    mutate(subject_to_remove = ifelse(wt2 > wt1, fid1_iid1, fid2_iid2))
  
  if(phenotypic_naive){
    singular_nodes = singular_nodes %>%
      mutate(subject_to_remove = ifelse(runif(1>0.5), fid1_iid1, fid2_iid2))
  }
  
  subject_to_remove_list = c(subject_to_remove_list, singular_nodes[["subject_to_remove"]])
  
  # Report
  kinship = kinship %>%
    filter(!(fid1_iid1 %in% subject_to_remove_list)) %>%
    filter(!(fid2_iid2 %in% subject_to_remove_list))
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  average_wt = (pheno %>%
                  filter(!(fid_iid %in% subject_to_remove_list)))[["wt"]] %>%
    mean()
  cat(paste0("Average subject phenotypic weight: ", round(average_wt,2), "\n\n"))
  
  n_remove = length(subject_to_remove_list)
  n_total = dim(pheno)[1]
  cat(paste0("Suggest removing ", n_remove, " subjects out of ", n_total, ".\n"))
  
  subject_to_remove_list = tibble(
    FID = strsplit(subject_to_remove_list, split = "_") %>%
      lapply(function(x){x[2]}) %>% unlist(),
    IID = strsplit(subject_to_remove_list, split = "_") %>%
      lapply(function(x){x[3]}) %>% unlist()
  )
  
  return(subject_to_remove_list)
  
}

### NOT RUN
# subject_to_remove = kdps(phenotype_file = "simulation/pheno.txt",
#                          kinship_file = "simulation/kinship.txt",
#                          fuzziness = 0,
#                          phenotype_name = "pheno2",
#                          prioritize_high = FALSE,
#                          prioritize_low = FALSE,
#                          phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY"),
#                          # Header names in the phenotype file
#                          fid_name = "FID",
#                          iid_name = "IID",
#                          # Header names in the kinship file
#                          fid1_name = "FID1",
#                          iid1_name = "IID1",
#                          fid2_name = "FID2",
#                          iid2_name = "IID2",
#                          kinship_name = "KINSHIP",
#                          kinship_threshold = 0.0442)

