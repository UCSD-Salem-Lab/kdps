library(dplyr)
library(tidyr)
library(data.table)
library(tibble)

phenotype_file = "simulation/pheno.txt"
fid_name = "FID"
iid_name = "IID"
phenotype_name = "pheno2"
prioritize_high = FALSE
prioritize_low = FALSE
phenotype_rank = c("DISEASED1", "DISEASED2", "HEALTHY")

kinship_file = "simulation/kinship.txt"
fid1_name = "FID1"
iid1_name = "IID1"
fid2_name = "FID2"
iid2_name = "IID2"
kinship_name = "KINSHIP"
kinship_threshold = 0.0442

### Make phenotype data frame
pheno = as.data.frame(fread(phenotype_file))

if(!any(prioritize_high, prioritize_low)){
  phenotype = pheno[[phenotype_name]]
  mapping = setNames(1:length(phenotype_rank), phenotype_rank)
  wt = mapping[phenotype]
  names(wt) = NULL
  wt[is.na(wt)] = 0
}else{
  if(prioritize_high){
    phenotype = as.numeric(pheno[[phenotype_name]])
    na_wt = min(phenotype, na.rm = TRUE) - 1
    wt = ifelse(is.na(phenotype), na_wt, phenotype)
    wt = (wt - mean(wt))/sd(wt)
  }
  if(prioritize_low){
    phenotype = as.numeric(pheno[[phenotype_name]])
    na_wt = max(phenotype, na.rm = TRUE) + 1
    wt = ifelse(is.na(phenotype), na_wt, phenotype)
    wt = (wt - mean(wt))/sd(wt)
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




# pivot_wider(kinship, names_from = fid2_iid2, values_from = related)

subjects = unique(c(
  kinship[["fid1_iid1"]],
  kinship[["fid2_iid2"]]
))

kinship_mat = matrix(
  data = NA, 
  nrow = length(subjects),
  ncol = length(subjects)
)

intersect(kinship[["fid1_iid1"]], kinship[["fid2_iid2"]])

for(i in 1:dim(kinship_mat)[1]){
  for(j in 1:dim(kinship_mat)[2]){
    print(round(i/dim(kinship_mat)[1],3))
  }
}

relationship = table(c(
  kinship[["fid1_iid1"]],
  kinship[["fid2_iid2"]]
))

relationship = tibble(
  subject = names(relationship),
  count = as.vector(relationship)
) %>%
  arrange(desc(count))

kinship_mat = tibble(
  FID1_IID1 = paste("subject", as.character(unlist(data[,1])), as.character(unlist(data[,2])), sep = "_"),
  FID2_IID2 = paste("subject", as.character(unlist(data[,3])), as.character(unlist(data[,4])), sep = "_"),
  kinship = data[,which(grepl(pattern = "kinship", names(data), ignore.case = TRUE))]
)

kinship_mat = pivot_wider(kinship_mat, names_from = FID1_IID1, values_from = kinship)

kinship_mat[, match(kinship_mat[["FID1_IID1"]], colnames(kinship_mat))]

as.data.frame(data_cor) %>%
  rownames_to_column(var="subject1") %>%
  gather(key="subject2", value="kinship", -1) %>% dim()

data = as.data.frame(matrix(rnorm(n^2), nrow = n))
names(data) = paste0("SUBJECT", 1:n)

data_cor = cor(data)
data_cor[abs(data_cor) > 0.1] = 1
data_cor[abs(data_cor) <= 0.1] = 0
rownames(data_cor) = colnames(data_cor)

subject_weight = data.frame(
  subject_name = paste0("SUBJECT", 1:n),
  weight = round(runif(100, min = 0, max = n))
)




sum_omit_na = function(x){
  sum(x, na.rm = TRUE)
}

dimension = dim(data_cor)[1]
data_cor[diag(TRUE, nrow = dimension, ncol = dimension)] = NA

data_cor_margin = addmargins(data_cor, FUN = list(sum = sum_omit_na), quiet = TRUE)
keep_list = names(which(data_cor_margin[,dimension + 1][1:dimension] == 0))
data_cor = data_cor[, !(colnames(data_cor) %in% keep_list)]
data_cor = data_cor[!(rownames(data_cor) %in% keep_list), ]
dimension = dim(data_cor)[1]

remove_list = c()
remaining_relationships = sum(data_cor, na.rm = TRUE)/2

while(remaining_relationships != 0){
  data_cor_margin = addmargins(data_cor, FUN = list(sum = sum_omit_na), quiet = TRUE)
  max_relationship_count = max(data_cor_margin[,dimension + 1][1:dimension])
  super_subjects = names(which(data_cor_margin[,dimension + 1][1:dimension] == max_relationship_count))
  super_subjects_weight = subject_weight[subject_weight[["subject_name"]] %in% super_subjects,]
  subject_to_remove = super_subjects_weight[["subject_name"]][which.min(super_subjects_weight[["weight"]])]
  remove_list = c(remove_list, subject_to_remove)
  data_cor[, colnames(data_cor) == subject_to_remove] = NA
  data_cor[rownames(data_cor) == subject_to_remove, ] = NA
  remaining_relationships = sum(data_cor, na.rm = TRUE)/2
  
  cat(paste0("Removing ", subject_to_remove, "...\n"))
  cat(paste0(remaining_relationships, " relationships remaining...\n"))
  cat("\n\n")
}

keep_list = subject_weight[["subject_name"]][!subject_weight[["subject_name"]] %in% remove_list]

