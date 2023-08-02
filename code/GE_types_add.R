library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")


data <- read_csv("./results/2013-2022_approvals_GE_by_source.csv")
local_approvals <- read_csv("./data/2013-2022_approvals.csv")

# Conditions for any GE support for approval (score>0, interactor or related disease has GE)
data_hasAnyGE <- data %>% filter(noTarget == FALSE) %>% # remove all non-human targets
    mutate(score_bool = score > 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE),
           hasAnyGE = interactionAssociation_bool|phenotypeAssociation_bool|score_bool) %>%
    group_by(datasourceId, originalDrugName, yearApproval) %>%
    summarise(hasAnyGE = any(hasAnyGE)) %>% 
    ungroup() %>%
    group_by(originalDrugName, yearApproval) %>%
    summarise(hasAnyGE = any(hasAnyGE)) %>% 
    ungroup() %>%
    group_by(yearApproval)

# GE support subtyping by evidence (original target/interactor + original disease/related disease)
data_GE_types <- data_hasAnyGE %>% 
    mutate(score_zero_bool = score == 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE)
           ) %>%
    group_by(Drug_name_original, yearApproval) %>%
    summarise(
      hasInteractorDiseasesGE = any(interactionAssociation_bool & score_zero_bool) & !any(!score_zero_bool),
      hasTargetRelatedGE = any(phenotypeAssociation_bool & score_zero_bool) & !any(!score_zero_bool),
      hasInteractorRelatvedGE = any(phenotypeAssociation_bool & interactionAssociation_bool & score_zero_bool) & !any(!score_zero_bool),
      hasTargetDiseaseGE = any(hasAnyGE) & !any(hasInteractorDiseasesGE) & !any(hasTargetRelatedGE) & !any(hasInteractorRelatvedGE)
    ) %>%
    ungroup()

approvals_GE <- local_approvals %>%
  # if complex drug has non-human target and human target, has_GE = NA for all of them
  left_join(data_GE_types, by = "originalDrugName") %>%
  mutate(hasAnyGE = ifelse(is.na(targetIds), NA, hasAnyGE)) %>%
  group_by(originalDrugName) %>%
  mutate(hasAnyGE = ifelse(is.na(max(hasAnyGE)), NA, hasAnyGE)) %>%
  ungroup()

write.table(approvals_GE, sep = ",", file = "./results/2013-2022_approvals_GE.csv", row.names = FALSE)