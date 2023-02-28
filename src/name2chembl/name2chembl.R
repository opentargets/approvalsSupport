library("tidyverse")
library("sparklyr")
library("sparklyr.nested")


data_release <- "22.11"

# Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint

# spark connect
sc <- spark_connect(master = "yarn", config = config)

# Approvals as reported in NRDD article (complex drugs should be listed using " ;")
local_next_approvals <- read_csv("./data/2018-2020/FDA_2018-2020_0.csv")
approvals_next_init <- sdf_copy_to(sc, local_next_approvals, overwrite = TRUE)

# Split and explode multiple Drug
approvals_next <- approvals_next_init %>%
    mutate(Drug = split(as.character(Drug), "; ")) %>%
    sdf_explode(Drug, keep_all = TRUE) 

#approvals_next %>% 
#    write.csv("Exploded_next.csv")

# Read Platform data (Dataset: Drug)
gs_path <- "gs://open-targets-data-releases/"
drug_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/molecule/",
    sep = ""
)

# ChEMBL ID from Platform
drug_info <- spark_read_parquet(sc, drug_path) %>%
    select(id, name) %>%
    sdf_distinct() 

#drug_info %>% 
#    write.csv("drug_info.csv")


# Drug name to chembl id
approvals_chembl <- approvals_next %>%
    mutate("Drug" = toupper(Drug)) %>% 
    left_join(drug_info, by = c("Drug" = "name")) %>% 
    collect()

approvals_chembl %>% 
    write.csv("name2chembl.csv")
