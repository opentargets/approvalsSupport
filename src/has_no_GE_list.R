library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

options(dplyr.width = Inf)
data <- read_csv("./output/2018-2022_approvals_v2.2.csv")
local_approvals <- read_csv("./data/2018-2022/2018-2022_approvals_v2.csv")

# Add MoAs to final table
# Extra MoAs required to fill the gaps
new_moas <- read_csv("./data/ammend_data/amend_moas.csv")
new_moas <- sdf_copy_to(sc, new_moas, overwrite = TRUE)

# available MoAs + ammended
moa <- spark_read_parquet(sc, moa_path, memory = FALSE) %>%
    select(chemblIds, targets) %>%
    sdf_explode(chemblIds) %>%
    sdf_explode(targets) %>%
    rename(targetId = targets) %>%
    sdf_distinct() %>%
    sdf_bind_rows(new_moas)

MoAs <- local_approvals %>%
  left_join(moa, by = c("DrugId" = "chemblIds"), copy=TRUE) %>%
  mutate(targetId = replace_na(targetId, "")) %>%
  group_by(Drug_name) %>%
  summarise(targetIds = paste(targetId, collapse = "; "))

local_approvals_MoAs <- local_approvals %>%
  left_join(MoAs, by = "Drug_name")

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

approvals_GE <- local_approvals_MoAs %>%
  left_join(data_humT_metadata, by = "Drug_name_original", copy = TRUE)

write.table(approvals_GE, sep = ",", file = "./output/2018-2022_approvals_GE_v2.2.csv", row.names = FALSE)


