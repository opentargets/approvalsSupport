library("tidyverse")
library("sparklyr")
library("sparklyr.nested")


data_release <- "23.02"

# Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint

# spark connect
sc <- spark_connect(master = "yarn", config = config)


local_approvals <- read_csv("./data/2013-2022/v3_full/2013-2022_approvals_GE_v3.2_in.csv")


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
  summarise(targetIds = paste(targetId, collapse = "; ")) %>%
  rowwise() %>%
  mutate(targetIds = toString(unique(unlist(strsplit(targetIds, "; ")))))

local_approvals_MoAs <- local_approvals %>%
  left_join(MoAs, by = "Drug_name") 


write.table(local_approvals_MoAs , sep = ",", file = "./data/2013-2022/v3_full/2013-2022_approvals_GE_v3.3_in.csv", row.names = FALSE)

