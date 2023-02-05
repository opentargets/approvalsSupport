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

# Approvals as reported in NRDD article
local_approvals <- read_csv("./data/all_approvals.csv")
approvals_init <- sdf_copy_to(sc, local_approvals, overwrite = TRUE)

# Split and explode multiple DrugId
approvals <- approvals_init %>%
    mutate( OriginalDrugId = DrugId,
            DrugId = split(as.character(DrugId), "; ")) %>%
    sdf_explode(DrugId, keep_all = TRUE) 

# Split and explode multiple DiseaseId
approvals <- approvals %>%
    mutate(DiseaseId = split(as.character(DiseaseId), "; ")) %>%
    sdf_explode(DiseaseId, keep_all = TRUE) 

# approvals %>% 
#     write.csv("Exploded.csv")

# Datasource metadata
local_ds_metadata <- read_csv("./data/datasourceMetadata.csv")
ds_names <- sdf_copy_to(sc, local_ds_metadata, overwrite = TRUE) %>% collect()

# Read Platform data
gs_path <- "gs://open-targets-data-releases/"
all_evidence_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/evidence/",
    sep = ""
)
moa_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/mechanismOfAction/",
    sep = ""
)
ass_indirectby_ds_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/associationByDatasourceIndirect/",
    sep = ""
)
disease_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/diseases/",
    sep = ""
)
interaction_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/interaction/",
    sep = ""
)
disease2phenotype_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/diseaseToPhenotype/",
    sep = ""
)

# Mechanisms of action
# Extra MoAs required to fill the gaps
new_moas <- read_csv("./data/amend_moas.csv")
new_moas <- sdf_copy_to(sc, new_moas, overwrite = TRUE)

# available MoAs + ammended
moa <- spark_read_parquet(sc, moa_path, memory = FALSE) %>%
    select(chemblIds, targets) %>%
    sdf_explode(chemblIds) %>%
    sdf_explode(targets) %>%
    rename(targetId = targets) %>%
    sdf_distinct() %>%
    sdf_bind_rows(new_moas)

# Platform ssociations indirect (by datasource)
ass_indirectby_ds <- spark_read_parquet(sc, ass_indirectby_ds_path)

# Joining associations information
ass <- approvals %>%
    rename(diseaseId = DiseaseId) %>%
    left_join(moa, by = c("DrugId" = "chemblIds")) %>%
    left_join(ass_indirectby_ds, by = c("diseaseId", "targetId")) %>%
    collect()

# Data about molecular interactions
interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
    filter(sourceDatabase == "intact") %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(scoring > 0.42) %>%
    select(targetA, targetB) %>%
    sdf_distinct()

interactors_ass <- approvals %>%
    rename(diseaseId = DiseaseId) %>%
    inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(interactions, by = c("targetId" = "targetA")) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("diseaseId" = "diseaseId", "targetB" = "targetId")
    ) %>%
    select(datasourceId, Drug_brand_name) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(interactionAssociation = TRUE)

# Additional phenotype curation
ammend_phenotypes <- read_csv("./data/amend_phenotypes.csv")
new_phenotypes <- sdf_copy_to(sc, ammend_phenotypes, overwrite = TRUE)

# Platform disease to phenotype data
disease2phenotype <- spark_read_parquet(
    sc,
    disease2phenotype_path,
    memory = FALSE
) %>%
    select(diseaseId = disease, phenotype) %>%
    sdf_distinct()

# Associations through indirect phenotypes
phenotype_ass <- approvals %>%
    rename(diseaseId = DiseaseId) %>%
    inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(
        disease2phenotype %>%
            sdf_bind_rows(new_phenotypes),
        by = c("diseaseId")
    ) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("phenotype" = "diseaseId", "targetId")
    ) %>%
    select(datasourceId, Drug_brand_name) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(phenotypeAssociation = TRUE)

# Data to plot
data2plot <- ass %>%
    select(datasourceId, Drug_brand_name, score) %>%
    complete(datasourceId, Drug_brand_name) %>%
    mutate(score = replace_na(score, 0)) %>%
    filter(!is.na(datasourceId)) %>%
    # TA
    left_join(
        ass %>%
            select(
                Drug_brand_name,
                TA
            ) %>%
            distinct(),
        by = "Drug_brand_name"
    ) %>%
    # targets
    left_join(
        ass %>%
            mutate(noTarget = is.na(targetId)) %>%
            select(
                Drug_brand_name,
                noTarget
            ) %>%
            distinct(),
        by = "Drug_brand_name"
    ) %>%
    # interactions
    left_join(
        interactors_ass,
        by = c("datasourceId", "Drug_brand_name")
    ) %>%
    # related phenotypes
    left_join(
        phenotype_ass,
        by = c("datasourceId", "Drug_brand_name")
    ) %>%
    mutate(
        interactionAssociation = ifelse(score > 0, TRUE, interactionAssociation)
    ) %>%
    mutate(
        phenotypeAssociation = ifelse(score > 0, TRUE, phenotypeAssociation)
    ) %>%
    mutate(score = ifelse(noTarget, NA, score)) %>%
    mutate(TA = ifelse(noTarget, "No human target", TA)) %>%
    mutate(
        TA = fct_other(
            TA,
            keep = c("Oncology", "No human target"),
            other_level = "Other indication"
        )
    ) %>%
    mutate(
        TA = fct_relevel(TA, c(
            "Oncology",
            "Other indication",
            "No human target"
        ))
    ) %>%
    # mutate(datasourceId = fct_relevel(datasourceId, names(ds_name_list))) %>%
    filter(!(datasourceId %in% c("chembl", "expression_atlas", "sysbio", "europepmc", "phenodigm", "reactome", "phewas_catalog"))) %>%
    # drug score for the purpose of reordering them
    mutate(rankscore = replace_na(score, 0)) %>%
    mutate(rankscore = ifelse(!is.na(interactionAssociation), rankscore + 0.01, rankscore)) %>%
    mutate(rankscore = ifelse(!is.na(phenotypeAssociation), rankscore + 0.03, rankscore)) %>%
    mutate(Drug_brand_name = fct_rev(fct_reorder(
        Drug_brand_name, rankscore, mean,
        na.rm = TRUE, .desc = TRUE
    ))) %>%
    group_by(
        datasourceId,
        Drug_brand_name,
        TA,
        noTarget,
        interactionAssociation,
        phenotypeAssociation
    ) %>%
    summarise(score = suppressWarnings(max(score, na.rm = TRUE))) %>%
        mutate(score = ifelse(score < 0, NA, score)) %>%
        left_join(ds_names, by = "datasourceId") %>%
        mutate(
            datasourceName = factor(datasourceName, levels = ds_names$datasourceName),
            datasourceType = factor(datasourceType, levels = c("Somatic", "Functional genomics (cancer)", "Rare mendelian", "Common disease"))
        ) %>%
        left_join(approvals %>% select(Drug_brand_name, Year) %>% collect(), by = "Drug_brand_name")

write.table(data2plot, sep = ",", file = "./output/approvals_support.csv", row.names = FALSE)
