library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

options(dplyr.width = Inf)
data <- read_csv("./output/v3/2013-2022_approvals_GE_v3.2_out.csv")
local_approvals <- read_csv("./output/v3/2013-2022_approvals_GE_v3.4_out.csv")


data_humT_metadata <- data %>% filter(noTarget == FALSE) %>% # remove all non-human targets
    mutate(score_bool = score > 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE),
           has_GE_ds = case_when(
            (TA == "Oncology" & interactionAssociation_bool|phenotypeAssociation_bool|score_bool) ~ TRUE,
            (datasourceType %in% c("Rare mendelian", "Common disease", "Mouse model") & TA %in% c("Other indication") 
            & interactionAssociation_bool|phenotypeAssociation_bool|score_bool) ~ TRUE,
            TRUE ~ FALSE
          )) %>%
    group_by(datasourceId, Drug_name_original, Year) %>%
    summarise(has_GE_ds = any(has_GE_ds)) %>% 
    ungroup() %>%
    group_by(Drug_name_original, Year) %>%
    summarise(has_GE_ds = any(has_GE_ds)) %>% 
    ungroup() %>%
    group_by(Year)

approvals_GE <- local_approvals %>%
  # if complex drug has non-human target and human target, has_GE = NA for all of them
  left_join(data_humT_metadata, by = "Drug_name_original") %>%
  mutate(has_GE_ds = ifelse(is.na(targetIds), NA, has_GE_ds)) %>%
  group_by(Drug_name_original) %>%
  mutate(has_GE_ds = ifelse(is.na(max(has_GE_ds)), NA, has_GE_ds)) %>%
  ungroup() %>%
  group_by(Brand_name, Drug_name_original, Sponsor, Properties, TA,
           Indication, Indication_EFO, DiseaseId, Year = Year.x, has_GE, has_GE_ds) %>%
  summarise(Drug_name = paste(Drug_name, collapse = "; "), 
            DrugId = paste(DrugId, collapse = "; "),
            targetIds = paste(targetIds, collapse = "; ")) %>%
  ungroup()

# approvals_GE <- aggregate(cbind(Drug_name, DrugId, targetIds) ~ Brand_name, data = approvals_GE, FUN = function(x) paste(x, collapse = "; "))

write.table(approvals_GE, sep = ",", file = "./output/v3/2013-2022_approvals_GE_v3.6_out.csv", row.names = FALSE)



