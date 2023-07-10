library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

options(dplyr.width = Inf)
data <- read_csv("./output/v3/2013-2022_approvals_GE_v3.2_out.csv")
local_approvals <- read_csv("./output/v3/2013-2022_approvals_GE_v6_novelty.csv")


data_humT_metadata <- data %>% filter(noTarget == FALSE) %>% # remove all non-human targets
    mutate(score_bool = score == 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE)
           ) %>%
    group_by(Drug_name_original, Year) %>%
    summarise(
      Interactor = any(interactionAssociation_bool & score_bool) & !any(!score_bool),
      Related = any(phenotypeAssociation_bool & score_bool) & !any(!score_bool),
      Interactor_Related = any(phenotypeAssociation_bool & interactionAssociation_bool & score_bool) & !any(!score_bool) 
    ) %>%
    ungroup()

approvals_GE <- local_approvals %>%
  # if complex drug has non-human target and human target, has_GE = NA for all of them
  left_join(data_humT_metadata, by = "Drug_name_original")

write.table(approvals_GE, sep = ",", file = "./output/v3/2013-2022_approvals_GE_v6.2_novelty.csv", row.names = FALSE)



new_data <- read_csv("./output/v3/2013-2022_approvals_GE_v6.3_novelty.csv")  
        # %>%
        # mutate(diseaseId = strsplit(diseaseId, "; "))
        # %>%
        # unnest(diseaseId)

amend_phenotypes <- read_csv("./data/ammend_data/amend_phenotypes_v2.csv") %>%
        group_by(diseaseId) %>%
        summarise(phenotype = paste(phenotype, collapse = ",")) %>%
        ungroup()

data_amend <- new_data %>% 
      left_join(amend_phenotypes, by = "diseaseId")

write.table(data_amend, sep = ",", file = "./output/v3/2013-2022_approvals_GE_v6.4_novelty.csv", row.names = FALSE)
