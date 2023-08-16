library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("dplyr")


data_release <- "23.02"

# Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint

# Spark connect
sc <- spark_connect(master = "yarn", config = config)

# Renaming to align with previous analysis (2021 approved drugs by D. Ochoa)
input_string <- "
inner_name,outer_name
Year,yearApproval
DrugId,drugId
Brand_name,brandDrugName
Drug_name_original,originalDrugName
Drug_name,drugName
targetId_original,targetIds
targetId_interactor,interactorIds
targetId_type,targetType
diseaseId,diseaseIds
diseaseId_related,relatedIds
Indication,indication
diseaseId_type,diseaseType
Disease_class,diseaseClass
has_GE,hasTargetDiseaseGE
Interactor,hasInteractorDiseaseGE
Related,hasTargetRelatedGE
Interactor_Related,hasInteractorRelatedGE
Sponsor,sponsor
Properties,properties
TA,therapeuticArea
Review_type,reviewType
"

data <- read.table(text = input_string, header = TRUE, stringsAsFactors = FALSE, sep = ",")
data_tibble <- as_tibble(data)
inside_rename_mapping <- setNames(data_tibble$outer_name, data_tibble$inner_name)
outside_rename_mapping <- setNames(data_tibble$inner_name, data_tibble$outer_name)


# Approvals as reported in NRDD article
local_approvals <- read_csv("./data/2013-2022_approvals_in.csv")
approvals_init <- sdf_copy_to(sc, local_approvals, overwrite = TRUE)

approvals_init <- approvals_init %>% 
    rename(any_of(c(!!!inside_rename_mapping)))

# Split and explode multiple DrugId
approvals <- approvals_init %>%
    mutate(OriginalDrugId = DrugId,
            DrugId = split(as.character(DrugId), ",")) %>%
    sdf_explode(DrugId, keep_all = TRUE) 

# Split and explode multiple DiseaseId
approvals <- approvals %>%
    mutate(diseaseId = split(as.character(diseaseId), ",")) %>%
    sdf_explode(diseaseId, keep_all = TRUE) 

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
    inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(interactions, by = c("targetId" = "targetA")) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("diseaseId" = "diseaseId", "targetB" = "targetId")
    ) %>%
    select(datasourceId, DrugId) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(interactionAssociation = TRUE)

# Additional phenotype curation
ammend_phenotypes <- read_csv("./data/amendPhenotypes.csv")
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
    inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(
        disease2phenotype %>%
            sdf_bind_rows(new_phenotypes),
        by = c("diseaseId")
    ) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("targetId", "phenotype" = "diseaseId") 
    ) %>%
    select(datasourceId, DrugId) %>%
    sdf_distinct() %>%
    collect() %>%
    mutate(phenotypeAssociation = TRUE)

# Data to plot
data2plot <- ass %>%
    select(datasourceId, DrugId, score) %>%
    complete(datasourceId, DrugId) %>%
    mutate(score = replace_na(score, 0)) %>%
    filter(!is.na(datasourceId)) %>%
    # TA
    left_join(
        ass %>%
            select(
                DrugId,
                TA
            ) %>%
            distinct(),
        by = "DrugId"
    ) %>%
    # targets
    left_join(
        ass %>%
            mutate(noTarget = is.na(targetId)) %>%
            select(
                DrugId,
                noTarget
            ) %>%
            distinct(),
        by = "DrugId"
    ) %>%
    # interactions
    left_join(
        interactors_ass,
        by = c("datasourceId", "DrugId")
    ) %>%
    # related phenotypes
    left_join(
        phenotype_ass,
        by = c("datasourceId", "DrugId")
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
    mutate(DrugId = fct_rev(fct_reorder(
        DrugId, rankscore, mean,
        na.rm = TRUE, .desc = TRUE
    ))) %>%
    group_by(
        datasourceId,
        DrugId,
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
            datasourceType = factor(datasourceType, levels = c("Somatic", "Functional genomics (cancer)", "Rare mendelian", "Common disease", "Mouse model"))
        ) %>%
        left_join(approvals %>% select(DrugId, Year, Brand_name,	Drug_name_original) %>% collect(), by = "DrugId")

 data2plot <- data2plot %>% 
    rename(any_of(c(!!!outside_rename_mapping)))

write.table(data2plot, sep = ",", file = "./results/2013-2022_approvals_GE_src.csv", row.names = FALSE)


# Add data about targets (targetIds)
approvals_moas <- approvals %>%
    left_join(moa %>% select(chemblIds, targetId), by = c("DrugId" = "chemblIds")) %>%
    rename(targetId_original = targetId)

# Add data about interactors (interactorIds)
approvals_inter <- approvals_moas %>%
    # inner_join(moa, by = c("drugId" = "chemblIds")) %>%
    left_join(interactions, by = c("targetId_original" = "targetA")) %>%
    left_join(ass_indirectby_ds %>% select(diseaseId, targetId), by = c("diseaseId", "targetB" = "targetId")) %>%
    rename(targetId_interactor = targetB)

# Add data about related conditions (relatedIds)
approvals_related <- approvals_inter %>%
  left_join(new_phenotypes %>% select(diseaseIds, phenotype), by = c("diseaseId" = "diseaseIds")) %>%
  rename(diseaseId_related = phenotype) %>%
  rename(any_of(c(!!!outside_rename_mapping))) %>%
  collect()

 approvals_final <- approvals_related %>% 
    group_by(brandDrugName) %>%
    summarise(across(everything(), ~ paste(unique(.x), collapse = ",")), .groups = "drop")

write.table(approvals_final, sep = ",", file = "./results/2013-2022_approvals_GE_out.csv", row.names = FALSE)
