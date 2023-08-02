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

# Initial testset of 2013-2022 approvals, exploded by originalDrugName with corresponding ChEMBL IDs
local_approvals <- read_csv("./data/2013-2022_approvals_input.csv")
approvals <- sdf_copy_to(sc, local_approvals, overwrite = TRUE)

# Split and explode multiple DrugId
approvals_drugs_expl <- approvals %>%
        mutate(drugId = split(as.character(drugId), ",")) %>%
        sdf_explode(drugId, keep_all = TRUE)

# Split and explode multiple diseaseIds
approvals_dis_expl <- approvals_drugs_expl %>%
    mutate(diseaseIds = split(as.character(diseaseIds), ",")) %>%
    sdf_explode(diseaseIds, keep_all = TRUE) 

# Datasource metadata
local_ds_metadata <- read_csv("./data/datasourceMetadata.csv")
ds_names <- sdf_copy_to(sc, local_ds_metadata, overwrite = TRUE) %>% collect()

# Read Platform data
gs_path <- "gs://open-targets-data-releases/"
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
new_moas <- read_csv("./data/amendMoas.csv")
new_moas <- sdf_copy_to(sc, new_moas, overwrite = TRUE)

# available + ammended MoAs
moa <- spark_read_parquet(sc, moa_path, memory = FALSE) %>%
    select(chemblIds, targets) %>%
    sdf_explode(chemblIds) %>%
    sdf_explode(targets) %>%
    rename(targetIds = targets) %>%
    sdf_distinct() %>%
    sdf_bind_rows(new_moas)

# Collect available + ammended MoAs for approvals
moas <- approvals_dis_expl %>%
  left_join(moa, by = c("drugId" = "chemblIds"), copy = TRUE) %>%
  group_by(drugName) %>%
  summarise(targetIds = concat_ws(",", collect_list(targetIds)))

# Add available + MoAs ammended to approvals dataset
approvals_moas <- approvals_dis_expl %>%
    left_join(moas, by = "drugName") 

# Platform associations indirect (by datasource)
ass_indirectby_ds <- spark_read_parquet(sc, ass_indirectby_ds_path)

# Joining associations information
ass <- approvals_moas %>%
    # rename(diseaseId = diseaseIds) %>%
    left_join(ass_indirectby_ds, by = c("diseaseIds" = "diseaseId", "targetIds" = "targetId")) %>%
    collect()

# Data about molecular interactions
interactions <- spark_read_parquet(sc, interaction_path, memory = FALSE) %>%
    filter(sourceDatabase == "intact") %>%
    filter(!is.na(targetA)) %>%
    filter(!is.na(targetB)) %>%
    filter(scoring > 0.42) %>%
    select(targetA, targetB) %>%
    sdf_distinct()

interactors_ass <- approvals_moas %>%
    # rename(diseaseId = DiseaseId) %>%
    # inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(interactions, by = c("targetIds" = "targetA")) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("diseaseIds" = "diseaseId", "targetB" = "targetId")
    ) %>%
    select(datasourceId, drugName, targetB) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(interactionAssociation = TRUE)


# Add data about molecular interactions to approvals dataset
approvals_inter <- approvals_moas %>%
    # inner_join(moa, by = c("drugId" = "chemblIds")) %>%
    inner_join(interactions, by = c("targetIds" = "targetA")) %>%
    inner_join(ass_indirectby_ds, by = c("diseaseIds" = "diseaseId", "targetB" = "targetId")) %>%
    rename(interactorIds = targetB)

# Additional phenotype curation
ammend_phenotypes <- read_csv("./data/amendPhenotypes.csv")
new_phenotypes <- sdf_copy_to(sc, ammend_phenotypes, overwrite = TRUE)

# Add data about related conditions to approvals dataset
approvals_related <- approvals_moas %>%
  left_join(new_phenotypes, by = "diseaseIds") %>%
  rename(relatedIds = phenotype)

# Platform disease to phenotype data
disease2phenotype <- spark_read_parquet(
    sc,
    disease2phenotype_path,
    memory = FALSE
) %>%
    select(diseaseIds = disease, phenotype) %>%
    sdf_distinct()

# Associations through indirect phenotypes
phenotype_ass <- approvals_related %>%
    # rename(diseaseId = diseaseIds) %>%
    # inner_join(moa, by = c("drugId" = "chemblIds")) %>%
    inner_join(
        disease2phenotype %>%
            sdf_bind_rows(new_phenotypes),
        by = c("diseaseIds")
    ) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("relatedIds" = "diseaseId", "targetIds" = "targetId")
    ) %>%
    select(datasourceId, drugName) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(phenotypeAssociation = TRUE)

# Data to plot
data2plot <- ass %>%
    select(datasourceId, drugName, score) %>%
    complete(datasourceId, drugName) %>%
    mutate(score = replace_na(score, 0)) %>%
    filter(!is.na(datasourceId)) %>%
    # therapeuticArea
    left_join(
        ass %>%
            select(
                drugName,
                therapeuticArea
            ) %>%
            distinct(),
        by = "drugName"
    ) %>%
    # targets
    left_join(
        ass %>%
            mutate(noTarget = is.na(targetIds)) %>%
            select(
                drugName,
                noTarget
            ) %>%
            distinct(),
        by = "drugName"
    ) %>%
    # interactions
    left_join(
        interactors_ass,
        by = c("datasourceId", "drugName")
    ) %>%
    # related phenotypes
    left_join(
        phenotype_ass,
        by = c("datasourceId", "drugName")
    ) %>%
    mutate(
        interactionAssociation = ifelse(score > 0, TRUE, interactionAssociation)
    ) %>%
    mutate(
        phenotypeAssociation = ifelse(score > 0, TRUE, phenotypeAssociation)
    ) %>%
    mutate(score = ifelse(noTarget, NA, score)) %>%
    mutate(therapeuticArea = ifelse(noTarget, "No human target", therapeuticArea)) %>%
    mutate(
        therapeuticArea = fct_other(
            therapeuticArea,
            keep = c("Oncology", "No human target"),
            other_level = "Other indication"
        )
    ) %>%
    mutate(
        therapeuticArea = fct_relevel(therapeuticArea, c(
            "Oncology",
            "Other indication",
            "No human target"
        ))
    ) %>%
    filter(!(datasourceId %in% c("chembl", "expression_atlas", "sysbio", "europepmc", "phenodigm", "reactome", "phewas_catalog"))) %>%
    # drug score for the purpose of reordering them
    mutate(rankscore = replace_na(score, 0)) %>%
    mutate(rankscore = ifelse(!is.na(interactionAssociation), rankscore + 0.01, rankscore)) %>%
    mutate(rankscore = ifelse(!is.na(phenotypeAssociation), rankscore + 0.03, rankscore)) %>%
    mutate(drugName = fct_rev(fct_reorder(
        drugName, rankscore, mean,
        na.rm = TRUE, .desc = TRUE
    ))) %>%
    group_by(
        datasourceId,
        drugName,
        therapeuticArea,
        noTarget,
        interactionAssociation,
        phenotypeAssociation
    ) %>%
    summarise(score = suppressWarnings(max(score, na.rm = TRUE))) %>%
        mutate(score = ifelse(score < 0, NA, score)) %>%
        left_join(ds_names, by = "datasourceId") %>%
        mutate(
            datasourceName = factor(datasourceName, levels = ds_names$datasourceName),
            datasourceType = factor(datasourceType, levels = c("Somatic", "Functional genomics (cancer)", "Rare mendelian", "Common disease", "Mouse model"))
        ) %>%
        left_join(approvals %>% select(drugName, yearApproval, brandDrugName, originalDrugName) %>% collect(), by = "drugName")

write.table(data2plot, sep = ",", file = "./results/2013-2022_approvals_GE_by_source_test.csv", row.names = FALSE)

# Prepare final approvals dataset
approvals_final <- approvals_inter %>%
  left_join(new_phenotypes, by = "diseaseIds") %>%
  rename(relatedIds = phenotype) %>%
  select(brandDrugName, originalDrugName, drugId, targetIds, interactorIds, sponsor, properties,
         therapeuticArea, indication, diseaseName, diseaseIds, diseaseClass, relatedIds,
         reviewType, yearApproval) %>%
  group_by(brandDrugName, originalDrugName, drugId) %>%
  summarise(
    drugId = concat_ws(",", drugId),
    targetIds = concat_ws(",", targetIds),
    interactorIds = concat_ws(",", interactorIds),
    diseaseIds = concat_ws(",", diseaseIds),
    relatedIds = concat_ws(",", relatedIds)
  )


write.table(approvals_final, sep = ",", file = "./results/2013-2022_approvals_GE_test.csv", row.names = FALSE)
