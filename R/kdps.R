#' Kinship Decouple and Phenotype Selection (KDPS)
#'
#' @description
#' The `kdps` function identifies subjects to be removed from a study based on kinship and phenotype
#' information. It uses kinship matrices and phenotype data to evaluate and prioritize subjects
#' according to their phenotype scores, taking into account their relatedness. The algorithm can
#' prioritize subjects with high or low phenotype values and filter out subjects based on kinship
#' thresholds and phenotype rankings. It aims to refine the study population by removing subjects
#' that do not meet specific genetic and phenotypic criteria, thus enhancing the robustness of
#' genetic association studies.
#'
#' @param phenotype_file A string specifying the path to the phenotype data file.
#' @param kinship_file A string specifying the path to the kinship matrix file.
#' @param fuzziness An integer representing the level of fuzziness allowed in removing related subjects, with a default of 0 (no fuzziness).
#' @param phenotype_name The name of the phenotype column in the phenotype file.
#' @param prioritize_high A logical indicating whether to prioritize subjects with high phenotype values for removal.
#' @param prioritize_low A logical indicating whether to prioritize subjects with low phenotype values for removal.
#' @param phenotype_rank A character vector specifying the ranking of phenotypes from highest priority (first) to lowest.
#' @param fid_name The column name for family IDs in the phenotype file.
#' @param iid_name The column name for individual IDs in the phenotype file.
#' @param fid1_name The column name for the first individual's family ID in the kinship file.
#' @param iid1_name The column name for the first individual's ID in the kinship file.
#' @param fid2_name The column name for the second individual's family ID in the kinship file.
#' @param iid2_name The column name for the second individual's ID in the kinship file.
#' @param kinship_name The name of the kinship score column in the kinship file.
#' @param kinship_threshold A numeric threshold for the kinship score, above which individuals are considered related.
#' @param phenotypic_naive A logical indicating whether to ignore phenotype information when resolving conflicts between related individuals.
#'
#' @details
#' The function first processes phenotype and kinship data from the specified files, then evaluates
#' subjects based on the provided parameters. It calculates weights for each subject based on their
#' phenotype and uses these weights along with the kinship information to identify subjects that
#' should be removed to minimize relatedness in the study population. The function offers flexibility
#' in handling phenotypes through ranking and prioritization options and can adjust the stringency of
#' relatedness filtering through the kinship threshold and fuzziness parameter.
#'
#' @return A data frame with two columns, `FID` and `IID`, representing the family and individual IDs
#' of subjects suggested for removal. This output can be used to refine the study population by
#' excluding these subjects in subsequent analyses.
#' @importFrom stats runif setNames
#' @examples
#' \dontrun{
#' subject_to_remove = kdps(phenotype_file = "data/simple_pheno.txt",
#'                          kinship_file = "data/simple_kinship.txt",
#'                          fuzziness = 0,
#'                          phenotype_name = "pheno2",
#'                          prioritize_high = FALSE,
#'                          prioritize_low = FALSE,
#'                          phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY"),
#'                          fid_name = "FID",
#'                          iid_name = "IID",
#'                          fid1_name = "FID1",
#'                          iid1_name = "IID1",
#'                          fid2_name = "FID2",
#'                          iid2_name = "IID2",
#'                          kinship_name = "KINSHIP",
#'                          kinship_threshold = 0.0442,
#'                          phenotypic_naive = FALSE)
#' }
#'
#' @export

kdps = function(phenotype_file = "data/simple_pheno.txt",
                kinship_file = "data/simple_kinship.txt",
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
  
  # Make phenotype data frame
  pheno = as.data.frame(data.table::fread(phenotype_file))
  pheno[[phenotype_name]] = as.character(pheno[[phenotype_name]])
  phenotype_rank = as.character(phenotype_rank)
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
  
  pheno = tibble::tibble(
    fid_iid = fid_iid,
    wt = wt
  )
  
  # Make kinship data frame
  kinship = as.data.frame(data.table::fread(kinship_file))
  kinship = tibble::tibble(
    fid1 = kinship[[fid1_name]],
    iid1 = kinship[[iid1_name]],
    fid2 = kinship[[fid2_name]],
    iid2 = kinship[[iid2_name]],
    kinship = kinship[[kinship_name]]
  ) 
  
  kinship = kinship |>
    dplyr::mutate(fid1_iid1 = paste(fid1, iid1, sep = "_")) |>
    dplyr::mutate(fid2_iid2 = paste(fid2, iid2, sep = "_")) |>
    dplyr::mutate(fid1_iid1 = paste0("subject_", fid1_iid1)) |>
    dplyr::mutate(fid2_iid2 = paste0("subject_", fid2_iid2)) |>
    dplyr::mutate(related = kinship >= kinship_threshold) |>
    dplyr::filter(related) |>
    dplyr::select(fid1_iid1, fid2_iid2, related)
  
  # Remove subjects that are not found in both the kinship and phenotype files
  kinship_subjects = unique(c(kinship[["fid1_iid1"]], kinship[["fid2_iid2"]]))
  pheno_subjects = unique(pheno[["fid_iid"]])
  common_subjects = intersect(kinship_subjects, pheno_subjects)
  
  cat(paste0(length(kinship_subjects), " unique related subjects found in the kinship file...\n"))
  cat(paste0(length(pheno_subjects), " unique subjects found in the phenotype file...\n"))
  cat(paste0(length(common_subjects), " subjects are found in both the related kinship and phenotype files...\n"))
  
  kinship_standout = kinship_subjects[!(kinship_subjects %in% common_subjects)]
  pheno_standout = pheno_subjects[!(pheno_subjects %in% common_subjects)]
  
  cat(paste0(length(kinship_standout), " unique subjects found only in the related kinship file...\n"))
  cat(paste0(length(pheno_standout), " unique subjects found only in the phenotype file...\n"))
  
  cat("Removing subjects from the kinship file without phenotype information...\n")
  subject_to_remove_list = kinship_standout
  
  kinship = kinship |>
    dplyr::filter(!(fid1_iid1 %in% subject_to_remove_list)) |>
    dplyr::filter(!(fid2_iid2 %in% subject_to_remove_list))
  
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  # Filter out singular relatedness
  cat("Filtering out singular relatedness...\n")
  relationship = table(c(
    kinship[["fid1_iid1"]],
    kinship[["fid2_iid2"]]
  ))
  one_timers = names(relationship[relationship == 1])
  
  singular_nodes = kinship |>
    dplyr::filter(fid1_iid1 %in% one_timers) |>
    dplyr::filter(fid2_iid2 %in% one_timers) |>
    dplyr::left_join(dplyr::mutate(pheno, fid1_iid1 = fid_iid, wt1 = wt) |> 
                       dplyr::select(fid1_iid1, wt1),
                     by = "fid1_iid1") |>
    dplyr::left_join(dplyr::mutate(pheno, fid2_iid2 = fid_iid, wt2 = wt) |> 
                       dplyr::select(fid2_iid2, wt2),
                     by = "fid2_iid2") |>
    dplyr::mutate(subject_to_remove = ifelse(wt2 > wt1, fid1_iid1, fid2_iid2))
  
  if(phenotypic_naive){
    singular_nodes = singular_nodes |>
      dplyr::mutate(subject_to_remove = ifelse(runif(dim(singular_nodes)[1]) > 0.5,
                                        fid1_iid1, fid2_iid2))
  }
  
  subject_to_remove_list = c(subject_to_remove_list, singular_nodes[["subject_to_remove"]])
  
  kinship = kinship |>
    dplyr::filter(!(fid1_iid1 %in% subject_to_remove_list)) |>
    dplyr::filter(!(fid2_iid2 %in% subject_to_remove_list))
  
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  if(dim(kinship)[1] == 0){
    cat("Only singular relatedness are present in the dataset...\n")
    average_wt = (pheno |>
                    dplyr::filter(!(fid_iid %in% subject_to_remove_list)))[["wt"]] |>
      mean()
    cat(paste0("Average subject phenotypic weight: ", round(average_wt,2), "\n\n"))
    
    subject_to_remove_list = tibble::tibble(
      FID = strsplit(subject_to_remove_list, split = "_") |>
        lapply(function(x){x[2]}) |> unlist(),
      IID = strsplit(subject_to_remove_list, split = "_") |>
        lapply(function(x){x[3]}) |> unlist()
    )
    
    return(subject_to_remove_list)
  }
  
  # Filter out singular nodes if fuzziness is 0
  if(fuzziness == 0){
    cat(paste0("Filtering out isolated super-subjects...\n"))
    relationship = table(c(
      kinship[["fid1_iid1"]],
      kinship[["fid2_iid2"]]
    ))
    isolated_n_timer = vector()
    for(n_relation in 2:10){
      cat(paste0("Filtering out isolated super-subjects with ", n_relation, " relatives...\n"))
      n_timers = names(relationship[relationship == n_relation])
      if(length(n_timers) > 0){
        pb = progress::progress_bar$new(total = length(n_timers))
        for(subject in n_timers){
          subject_node = kinship |>
            dplyr::filter(fid1_iid1 == subject | fid2_iid2 == subject)
          connected_subjects = c(subject_node[["fid1_iid1"]],
                                 subject_node[["fid2_iid2"]])
          connected_subjects = connected_subjects[connected_subjects != subject]
          if(all(connected_subjects %in% one_timers)){
            isolated_n_timer = rbind.data.frame(
              isolated_n_timer,
              tibble::tibble(subject = subject, n = n_relation)
            )
          }
          pb$tick()
        }
      }
    }
    
    subject_to_remove_list = c(subject_to_remove_list, 
                               isolated_n_timer[["subject"]])
    kinship = kinship |>
      dplyr::filter(!(fid1_iid1 %in% subject_to_remove_list)) |>
      dplyr::filter(!(fid2_iid2 %in% subject_to_remove_list))
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
    
    subject_to_remove = (tibble::tibble(
      fid_iid = names(relationship[relationship >= max_count_corrected])
    ) |>
      dplyr::left_join(pheno, by = "fid_iid") |>
      dplyr::arrange(wt))[["fid_iid"]][1]
    
    if(phenotypic_naive){
      subject_to_remove = (tibble::tibble(
        fid_iid = names(relationship[relationship >= max_count_corrected])
      ) |>
        dplyr::left_join(pheno, by = "fid_iid"))[["fid_iid"]]
      subject_to_remove = subject_to_remove[sample(1:length(subject_to_remove), 1)]
    }
    
    kinship = kinship |>
      dplyr::filter(fid1_iid1 != subject_to_remove) |>
      dplyr::filter(fid2_iid2 != subject_to_remove)
    
    cat(paste0("# of relationships left: ", dim(kinship)[1], "\r"))
    subject_to_remove_list = c(subject_to_remove_list, subject_to_remove)
  }
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  # Work up
  cat("Performing work-up...\n")
  singular_nodes = kinship |>
    dplyr::left_join(dplyr::mutate(pheno, fid1_iid1 = fid_iid, wt1 = wt) |> 
                       dplyr::select(fid1_iid1, wt1),
                     by = "fid1_iid1") |>
    dplyr::left_join(dplyr::mutate(pheno, fid2_iid2 = fid_iid, wt2 = wt) |> 
                       dplyr::select(fid2_iid2, wt2),
                     by = "fid2_iid2") |>
    dplyr::mutate(subject_to_remove = ifelse(wt2 > wt1, fid1_iid1, fid2_iid2))
  
  if(phenotypic_naive){
    singular_nodes = singular_nodes |>
      dplyr::mutate(subject_to_remove = ifelse(runif(1>0.5), fid1_iid1, fid2_iid2))
  }
  
  subject_to_remove_list = c(subject_to_remove_list, singular_nodes[["subject_to_remove"]])
  
  # Report
  kinship = kinship |>
    dplyr::filter(!(fid1_iid1 %in% subject_to_remove_list)) |>
    dplyr::filter(!(fid2_iid2 %in% subject_to_remove_list))
  cat(paste0("# of relationships left: ", dim(kinship)[1], "\n"))
  
  average_wt = (pheno |>
                  dplyr::filter(!(fid_iid %in% subject_to_remove_list)))[["wt"]] |>
    mean()
  cat(paste0("Average subject phenotypic weight: ", round(average_wt,2), "\n\n"))
  
  n_remove = length(subject_to_remove_list)
  n_total = dim(pheno)[1]
  cat(paste0("Suggest removing ", n_remove, " subjects out of ", n_total, ".\n"))
  
  subject_to_remove_list = tibble::tibble(
    FID = strsplit(subject_to_remove_list, split = "_") |>
      lapply(function(x){x[2]}) |> unlist(),
    IID = strsplit(subject_to_remove_list, split = "_") |>
      lapply(function(x){x[3]}) |> unlist()
  )
  
  return(subject_to_remove_list)
  
}

### NOT RUN
# subject_to_remove = kdps(phenotype_file = "data/simple_pheno.txt",
#                          kinship_file = "data/simple_kinship.txt",
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
#                          kinship_threshold = 0.0442,
#                          phenotypic_naive = FALSE)

