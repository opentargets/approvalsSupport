library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

options(dplyr.width = Inf)
data <- read_csv("./output/v3/2013-2022_approvals_GE_v3.2_out.csv")
local_approvals <- read_csv("./data/2013-2022/v3/2013-2022_approvals_GE_v3.2_in.csv")


data_humT_metadata <- data %>% filter(noTarget == FALSE) %>% # remove all non-human targets
    mutate(score_bool = score > 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE),
           has_GE = interactionAssociation_bool|phenotypeAssociation_bool|score_bool) %>%
    group_by(datasourceId, Drug_name_original, Year) %>%
    summarise(has_GE = any(has_GE)) %>% 
    ungroup() %>%
    group_by(Drug_name_original, Year) %>%
    summarise(has_GE = any(has_GE)) %>% 
    ungroup() %>%
    group_by(Year)

approvals_GE <- local_approvals %>%
  # if complex drug has non-human target and human target, has_GE = NA for all of them
  left_join(data_humT_metadata, by = "Drug_name_original") %>%
  mutate(has_GE = ifelse(is.na(targetIds), NA, has_GE)) %>%
  group_by(Drug_name_original) %>%
  mutate(has_GE = ifelse(is.na(max(has_GE)), NA, has_GE)) %>%
  ungroup()


write.table(approvals_GE, sep = ",", file = "./output/v3/2013-2022_approvals_GE_v3.3_out.csv", row.names = FALSE)



