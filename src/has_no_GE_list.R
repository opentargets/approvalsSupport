library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

options(dplyr.width = Inf)
data <- read_csv("./output/2013-2022_approvals_v1.2_out.csv")
local_approvals <- read_csv("./data/2013-2022/2013-2022_approvals_v1.2_trgts.csv")


data_humT_metadata <- data %>% filter(noTarget == FALSE) %>%
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
  # group_by(Drug_name_original, Year) %>%
  left_join(data_humT_metadata, by = "Drug_name_original") 
  # ungroup()

write.table(approvals_GE, sep = ",", file = "./output/2013-2022_approvals_GE_v1.2.csv", row.names = FALSE)


