library(dplyr)
library(tidyr)
library(data.table)
library(tibble)

set.seed(492357816)

n = 1000

data = as.data.frame(fread("test_kinship.txt"))
kinship_mat = tibble(
  FID1_IID1 = paste("SUBJECT", as.character(unlist(data[,1])), as.character(unlist(data[,2])), sep = "_"),
  FID2_IID2 = paste("SUBJECT", as.character(unlist(data[,3])), as.character(unlist(data[,4])), sep = "_"),
  kinship = data[,which(grepl(pattern = "kinship", names(data), ignore.case = TRUE))]
)

kinship_mat = pivot_wider(kinship_mat, names_from = FID2_IID2, values_from = kinship)

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

