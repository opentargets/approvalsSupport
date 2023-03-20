
 
directSources <- ass %>%
    filter(!(datasourceId %in% c("chembl", "expression_atlas", "sysbio", "europepmc", "phenodigm", "reactome", "phewas_catalog"))) %>%
    mutate(datasourceId = datasourceId %>% str_replace("eva", "clinvar")) %>%
    filter(!is.na(datasourceId)) %>%
    group_by(Drug_name) %>%
    summarise(directSources = paste(unique(datasourceId), collapse = ";"))

summaryResults <- output %>% 
    filter(datasourceType == "Any") %>%
    select(Drug_name, evidence)

closePhenotypes <- phenotype_ass %>% 
    select(Drug_name, datasourceId, phenotype) %>%
    left_join(
        spark_read_parquet(sc, disease_path) %>% 
        select(phenotype = id, phenotypeName = name),
        by = "phenotype") %>%
    collect() %>% 
    mutate(datasourceId = datasourceId %>% str_replace("eva", "clinvar")) %>%
    filter(!(datasourceId %in% c("chembl", "expression_atlas", "sysbio", "europepmc", "phenodigm", "reactome", "phewas_catalog"))) %>%
    distinct() %>% 
    group_by(Drug_name) %>%
    summarise(
        closePhenotypeIds = paste(unique(phenotype), collapse = ";"),
        closePhenotypeNames = paste(unique(phenotypeName), collapse = ";"),
        closePhenotypeDataSources = paste(unique(datasourceId), collapse = ";")
    )

target_path <- paste(
    gs_path, data_release,
    "/output/etl/parquet/target/",
    sep = ""
)

intDf <- approvals %>%
    rename(diseaseId = DiseaseId) %>%
    inner_join(moa, by = c("DrugId" = "chemblIds")) %>%
    inner_join(interactions, by = c("targetId" = "targetA")) %>%
    inner_join(
        ass_indirectby_ds,
        by = c("diseaseId" = "diseaseId", "targetB" = "targetId")
    ) %>%
    left_join(
        spark_read_parquet(sc, target_path) %>% 
        select(targetB = id, approvedSymbol),
        by = "targetB"
    ) %>%
    select(Drug_name, targetB, datasourceId, approvedSymbol) %>%
    collect() %>%
    mutate(datasourceId = datasourceId %>% str_replace("eva", "clinvar")) %>%
    filter(!(datasourceId %in% c("chembl", "expression_atlas", "sysbio", "europepmc", "phenodigm", "reactome", "phewas_catalog"))) %>%
    distinct() %>% 
    group_by(Drug_name) %>%
    summarise(
        interactingIds = paste(unique(targetB), collapse = ";"),
        interactingSymbols = paste(unique(approvedSymbol), collapse = ";"),
        interactingDataSources = paste(unique(datasourceId), collapse = ";")
    )

out <- ass %>% 
    group_by(Drug_name, Sponsor, DrugId, Indication, diseaseId, Properties) %>% 
    summarise(targetIds = paste(targetId, collapse = ";")) %>%
    left_join(summaryResults, by = "Drug_name") %>%
    left_join(directSources, by = "Drug_name") %>%
    left_join(closePhenotypes, by = "Drug_name") %>%
    left_join(intDf, by = "Drug_name")

out %>% write_csv("/home/ochoa/2021_approvals_output.csv")

