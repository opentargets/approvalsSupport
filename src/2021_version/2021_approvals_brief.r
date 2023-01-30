library("tidyverse")
library("sparklyr")
library("sparklyr.nested")
library("cowplot")
library("ggsci")

#Spark config
config <- spark_config()

# Allowing to GCP datasets access
config$spark.hadoop.fs.gs.requester.pays.mode <- "AUTO" # nolint
config$spark.hadoop.fs.gs.requester.pays.project.id <- "open-targets-eu-dev" # nolint

# spark connect
sc <- spark_connect(master = "local", config = config)

# Approvals as reported in NRDD article
gs_approvals <- "gs://ot-team/dochoa/2021_approvals.csv"
approvals <- spark_read_csv(
    sc,
    path = gs_approvals,
    memory = FALSE
)

# Datasource metadata
ds_names <- spark_read_csv(
    sc,
    path = "gs://ot-team/dochoa/datasourceMetadata.csv",
    memory = FALSE) %>%
    collect()

# Read Platform data
gs_path <- "gs://open-targets-data-releases/"
data_release <- "21.11"
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
ammend_moas <- list(
    "CHEMBL4594302" = "ENSG00000134318",
    "CHEMBL4297741" = "ENSG00000215644",
    "CHEMBL4297774" = "ENSG00000146648",
    "CHEMBL4297774" = "ENSG00000105976",
    "CHEMBL4298185" = "ENSG00000112964", # chembl missing in platform
    "CHEMBL4650319" = "ENSG00000146648",
    "CHEMBL1863514" = "ENSG00000166183",
    "CHEMBL4594320" = "ENSG00000171298"
)
new_moas <- data.frame(
    chemblIds = names(ammend_moas),
    targetId = unlist(ammend_moas)
)
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
ammend_phenotypes <- list(
    # Microalbuminuria (biomarker of CKD)
    "EFO_0000401" = "HP_0012594",
    # glycodeoxycholate sulfate (one of the bile acids that cause pruritus)
    "Orphanet_172" = "EFO_0005653",
    "Orphanet_52" = "EFO_0005653",
    # achondroplasia -> body height
    "Orphanet_15" = "EFO_0004339",
    "Orphanet_15" = "Orphanet_329191",
    #von hippel lindau -> renal carcinoma
    "Orphanet_892" = "EFO_0000681",
    "EFO_0001360" = "MONDO_0018582",
    # growth delay -> height
    "HP_0001510" = "EFO_0004339",
    #CAD -> myocardial infarctation
    "EFO_0001645" = "EFO_0000612"
)
new_phenotypes <- data.frame(
    diseaseId = names(ammend_phenotypes),
    phenotype = unlist(ammend_phenotypes)
)
new_phenotypes <- sdf_copy_to(sc, new_phenotypes, overwrite = TRUE)

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
        by = c("diseaseId")) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("phenotype" = "diseaseId", "targetId")) %>%
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
    #drug score for the purpose of reordering them
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
    )

# Values per data source
briefplotdata <- data2plot %>%
mutate(score = replace_na(score, 0)) %>%
group_by(Drug_brand_name, TA, datasourceType) %>%
summarise(
    noTarget = any(noTarget),
    interactionAssociation = any(interactionAssociation),
    phenotypeAssociation = any(phenotypeAssociation),
    score = ifelse(max(score, na.rm = TRUE) > 0, TRUE, FALSE)
) %>%
mutate(noTarget = replace_na(noTarget, FALSE)) %>%
mutate(phenotypeAssociation = replace_na(phenotypeAssociation, FALSE)) %>%
mutate(phenotypeAssociation = ifelse(score, FALSE, phenotypeAssociation)) %>%
mutate(interactionAssociation = replace_na(interactionAssociation, FALSE)) %>%
mutate(interactionAssociation = ifelse(score, FALSE, interactionAssociation)) %>%
mutate(interactionAssociation = ifelse(phenotypeAssociation, FALSE, interactionAssociation)) %>%
mutate(noEvidence = !(interactionAssociation | phenotypeAssociation | score | noTarget)) %>%
gather("evidence", "value", -Drug_brand_name, -TA, -datasourceType) %>%
filter(value)

# Values any data source
briefplotdataAny <- data2plot %>%
mutate(score = replace_na(score, 0)) %>%
group_by(Drug_brand_name, TA) %>%
summarise(
    noTarget = any(noTarget),
    interactionAssociation = any(interactionAssociation),
    phenotypeAssociation = any(phenotypeAssociation),
    score = ifelse(max(score, na.rm = TRUE) > 0, TRUE, FALSE)
) %>%
mutate(datasourceType = "Any") %>%
mutate(noTarget = replace_na(noTarget, FALSE)) %>%
mutate(phenotypeAssociation = replace_na(phenotypeAssociation, FALSE)) %>%
mutate(phenotypeAssociation = ifelse(score, FALSE, phenotypeAssociation)) %>%
mutate(interactionAssociation = replace_na(interactionAssociation, FALSE)) %>%
mutate(interactionAssociation = ifelse(score, FALSE, interactionAssociation)) %>%
mutate(interactionAssociation = ifelse(phenotypeAssociation, FALSE, interactionAssociation)) %>%
mutate(noEvidence = !(interactionAssociation | phenotypeAssociation | score | noTarget)) %>%
gather("evidence", "value", -Drug_brand_name, -TA, -datasourceType) %>%
filter(value)


output <- bind_rows(briefplotdataAny, briefplotdata) %>%
mutate(datasourceType = fct_relevel(datasourceType, levels = c("Any", "Somatic", "Functional genomics (cancer)", "Rare mendelian", "Common disease"))) %>%
mutate(evidence = fct_relevel(evidence,
    "score",
    "phenotypeAssociation",
    "interactionAssociation",
    "noTarget",
    "noEvidence")) %>%
mutate(evidence = fct_recode(evidence,
    "Direct" = "score",
    "Close phenotype" = "phenotypeAssociation",
    "Interacting protein" = "interactionAssociation",
    "No human target" = "noTarget",
    "Not available" = "noEvidence"
)) %>%
arrange(TA, desc(evidence)) %>%
group_by(datasourceType) %>% 
mutate(rn = row_number()) %>%
mutate(evidence = replace(evidence, evidence == "Not available", NA)) %>%
ggplot(aes(x = rn, y = fct_rev(datasourceType), fill = fct_rev(evidence))) +
geom_tile(color = "white", height = .8, size = 0.5) +
facet_grid(
    ~TA,
    scales = "free",
    space = "free"
) +
# scale_fill_npg(name = "Genetic support", na.value = "grey90") +
scale_fill_manual(
    name = "Genetic support",
    values = c("#3C5488FF", "#00A087FF", "#4DBBD5FF", "grey60"),
    breaks = c("Direct", "Close phenotype", "Interacting protein", "No human target"),
    na.value = "grey90") +
scale_y_discrete(name = "Genetic data source", labels = function(x) str_wrap(x, width = 12)) +
theme_cowplot(font_size = 11) +
theme(
    plot.background = element_rect(fill = "white"),
    strip.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.line = element_blank(),
    text = element_text(family = "sans"),
    panel.spacing = unit(-0.5, "lines")
)

ggsave(
    "/home/ochoa/2021_approvals_brief.pdf",
    plot = output,
    width = 6.5,
    height = 3.5,
    dpi = 400,
)
