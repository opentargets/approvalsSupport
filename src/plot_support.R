library("cowplot")
library("ggsci")
library(forcats)


data2plot <- read_csv("./output/2016-2017_approvals_v6.csv")


# Defining a function for plotting with parameters and one output by Year
ass_plotting <- function(year, version, width = 8, height = 11) {

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
            Drug_name,
            TA,
            interactionAssociation,
            phenotypeAssociation,
            Brand_name,
            Drug_name_original
        ) %>%
        gather("overlay", "value", -datasourceName, -datasourceType, -Drug_name, -TA, -Brand_name, -Drug_name_original) %>%
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
                fct_relevel(datasourceType, "Somatic", "Common disease", "Rare mendelian", "Functional genomics (cancer)", "Mouse model"), 
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
            legend.position = "bottom", #c(0.1, -0.4),
            legend.justification = c(0, 0),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank(),
            axis.line = element_blank(),
            text = element_text(family = "sans")
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

    print(paste0("filename is: ","./output/", year, "_approvals_", version, ".pdf"))

    ggsave(
        paste0("./output/", year, "_approvals_", version, ".pdf"),
        plot = output,
        width = width,
        height = height
    )
}

unique_years <- unique(data2plot$Year)

# Loop through unique years and apply function to relevant subset of data
for (year in unique_years) {
  ass_plotting(year, version = "v6") # apply function to subset
}
