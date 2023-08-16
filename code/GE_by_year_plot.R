library("cowplot")
library("ggsci")
library(forcats)
library("dplyr")
library("tidyverse")


data2plot <- read_csv("./results/2013-2022_approvals_GE_src.csv")

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

data2plot <- data2plot %>% 
    rename(any_of(c(!!!inside_rename_mapping)))

# Defining a function for plotting with parameters and one output by Year
ass_plotting <- function(year, version, width = 11, height = 12) {

data2plot$datasourceType <- as.character(data2plot$datasourceType) 
data2plot$datasourceType[data2plot$datasourceType == "Functional genomics (cancer)"] <- "Funct. gen.\n(cancer)"
data2plot$datasourceType[data2plot$datasourceType == "Common disease"] <- "Common\ndisease"
data2plot$datasourceType[data2plot$datasourceType == "Mouse model"] <- "Mouse\nmodel"
data2plot$datasourceType <- as.factor(data2plot$datasourceType)

# Complex drugs procession (if one of drug's compounds => non-human target, drug => non-human target)
    data2plot_proc <- data2plot %>%
        filter(Year == year) %>%
        mutate(has_no_human_target = TA == "No human target") %>%
        group_by(Drug_name_original) %>%
        summarise(guilty = max(has_no_human_target) - min(has_no_human_target)) %>%
        ungroup() %>%
        left_join(data2plot, by = "Drug_name_original") %>%
        mutate(TA = ifelse(guilty, "No human target", TA))

    # symbols to overlay in the plot
    overlay_data <- data2plot_proc %>%
        ungroup() %>%
        select(
            datasourceName,
            datasourceType,
            TA,
            interactionAssociation,
            phenotypeAssociation,
            Brand_name,
            Drug_name_original
        ) %>%
        gather("overlay", "value", -datasourceName, -datasourceType, -TA, -Brand_name, -Drug_name_original) %>%
        filter(!is.na(value)) %>%
        mutate(overlay = str_replace_all(overlay, "Association", "")) %>%
        mutate(overlaySize = ifelse(overlay == "phenotype", 3, 1)) %>%
        mutate(overlaySymbol = as.character(ifelse(overlay == "phenotype", 1, 16)))


    sort_order <- data2plot_proc %>%
    left_join(overlay_data) %>% 
    mutate(overlay_score = ifelse(value, 0.01, 0)) %>%
    group_by(Drug_name_original) %>%
    summarise(order = sum(score, na.rm=TRUE) + sum(overlay_score, na.rm=TRUE)) %>% arrange(-order) %>% ungroup()

    # plotting
    output <- data2plot_proc %>%
        left_join(sort_order) %>%
        mutate(Drug_name_original = fct_reorder(Drug_name_original, order)) %>%
        ggplot(aes(
            x = datasourceName,
            y = Drug_name_original
        )) +
        geom_tile(aes(fill = score), color = "white") +
        geom_point(
            data = overlay_data,
            aes(shape = overlay, size = overlaySize)
        ) +
        scale_fill_material("blue",
            na.value = "grey90",
            name = "Direct association"
        ) +
        scale_shape_manual(
            breaks = c("phenotype", "interaction"),
            labels = c("Direct or related phenotype", "Direct or interacting protein"),
            values = c(1, 16),
            name = "Supported by:"
        ) +
        scale_size_identity() +
        facet_grid(fct_relevel(TA, "Oncology", "Other indication", "No human target") ~
                fct_relevel(datasourceType, "Somatic", "Common\ndisease", "Rare mendelian", "Funct. gen.\n(cancer)", "Mouse\nmodel"), 
                scales = "free", space = "free") +
        theme_cowplot(font_size = 12) +
        # labs(
        #     title = "Supporting evidence on 2021 FDA drug approvals",
        #     subtitle = "Target-Disease evidence from Open Targets"
        #     # caption =
        #     #     "Source: Nat Reviews Drug Discovery 10.1038/d41573-022-00001-9"
        # ) +
        theme(
            plot.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            legend.direction = "horizontal",
            legend.box = "vertical",
            legend.position = c(-0.4, -0.15), #c(-0.45, -0.16), c(-0.5, -0.16)
            legend.justification = c(0, 0),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank(),
            axis.line = element_blank(),
            text = element_text(family = "sans"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 2.5), "cm")  # unit(c(0.5, 0.5, 0.5, 3), "cm"), unit(c(0.5, 0.5, 0.5, 5), "cm")
        ) +
        guides(
            fill = guide_colourbar(
                title.position = "top",
                title.hjust = 0.5,
                barwidth = 8,
                frame.colour = "black",
                ticks.colour = "black",
                order = 2
            ),
            shape = guide_legend(
                title.position = "top",
                direction = "vertical",
                order = 1
            )
        )

    print(paste0("filename is: ","./results/GE_by_year/", year, "_approvals_GE_", version, ".png"))

    ggsave(
        paste0("./results/GE_by_year/", year, "_approvals_GE_", version, ".png"),
        plot = output,
        width = width,
        height = height
    )
}

unique_years <- unique(data2plot$Year)

# Loop through unique years and apply function to relevant subset of data
for (year in unique_years) {
  ass_plotting(year, version = "5") # apply function to subset
}
