library("tidyverse")
library("sparklyr")
library("sparklyr.nested")


data_release <- "23.02"

# Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint


sc <- spark_connect(master = "yarn", config = config)


local_approvals <- read_csv("./output/2013-2022_approvals_v2_review.csv")
approvals <- sdf_copy_to(sc, local_approvals, overwrite = TRUE)


# Read Platform data (Dataset: Drug)
gs_path <- "gs://open-targets-data-releases/"
drug_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/molecule/",
    sep = ""
)

# ChEMBL ID from Platform
drug_info <- spark_read_parquet(sc, drug_path) %>%
    select(id, drugType) %>%
    sdf_distinct()


# Drug name to chembl id
approvals_mod <- approvals %>%
    left_join(drug_info, by = c("DrugId" = "id")) %>% 
    collect()

approvals_mod %>% 
    write.csv("./output/2013-2022_approvals_v2_review_mod.csv")