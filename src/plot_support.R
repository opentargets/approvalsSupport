library("cowplot")
library("ggsci")


data2plot <- read_csv("./output/approvals_support.csv")

data2plot <- data2plot %>% filter(Year == 2022)

# symbols to overlay in the plot
overlay_data <- data2plot %>%
    ungroup() %>%
    select(
        datasourceName,
        datasourceType,
        Drug_brand_name,
        TA,
        interactionAssociation,
        phenotypeAssociation
    ) %>%
    gather("overlay", "value", -datasourceName, -datasourceType, -Drug_brand_name, -TA) %>%
    filter(!is.na(value)) %>%
    mutate(overlay = str_replace_all(overlay, "Association", "")) %>%
    mutate(overlaySize = ifelse(overlay == "phenotype", 3, 1)) %>%
    mutate(overlaySymbol = as.character(ifelse(overlay == "phenotype", 1, 16)))

# plotting
output <- data2plot %>%
    ggplot(aes(
        x = datasourceName,
        y = Drug_brand_name
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
    facet_grid(TA ~ datasourceType, scales = "free", space = "free") +
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
        legend.position = c(-0.7, -0.16),
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

ggsave(
    "./output/2022_approvals_new.pdf",
    plot = output,
    width = 9,
    height = 15
)
