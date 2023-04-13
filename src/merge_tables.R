library(dplyr)
library("tidyverse")
library("cowplot")
library("ggsci")
library(forcats)
library("sparklyr")
library("sparklyr.nested")

data1 <- read_csv("./output/2016-2022_approvals_v1.csv")
data2 <- read_csv("./output/2014-2015_approvals_v2.2.csv")

local_approvals <- bind_rows(data1, data2)

write.table(local_approvals, sep = ",", file = "./output/2014-2022_approvals_v1.csv", row.names = FALSE)


data_release <- "23.02"

# Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint

# spark connect
sc <- spark_connect(master = "yarn", config = config)


# Add MoAs to final table
# Extra MoAs required to fill the gaps
new_moas <- read_csv("./data/ammend_data/amend_moas.csv")
new_moas <- sdf_copy_to(sc, new_moas, overwrite = TRUE)

gs_path <- "gs://open-targets-data-releases/"
moa_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/mechanismOfAction/",
    sep = ""
)

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
  left_join(MoAs, by = "Drug_name") %>%
  select(Brand_name,Drug_name_original,Drug_name,
        Sponsor,DrugId,Properties,
        Indication,Indication_EFO,DiseaseId,TA,
        Review_type,Year,targetIds = targetIds.y,
        has_GE,Why_no_GE,Related_indication_MoA,Related_DiseaseId_targetIds)
  

write.table(local_approvals_MoAs, sep = ",", file = "./data/2016-2022/2016-2022_approvals_v1.csv", row.names = FALSE)
