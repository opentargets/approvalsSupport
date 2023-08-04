library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(epitools)
library(ggpubr)

font = "Helvetica"

data <- read_csv("./results/2013-2022_approvals_GE_prec_cur.csv")

# Calculation of overall approval numbers per year
data_all_count <- data %>%
    select(originalDrugName, yearApproval) %>%
    distinct() %>%
    group_by(yearApproval) %>%
    summarise(n = n()) %>%
    ungroup()

# Creating a table with % of has_GE, has_no_GE and non-human targets
data_humT_metadata <- data %>% filter(hasAnyGE != "NA") %>%
        mutate(minDate = pmin(yearTargetDiseaseGE, yearInteractorDiseaseGE, yearTargetRelatedGE, yearInteractorRelatedGE, na.rm = TRUE),
                priorDate = ifelse(!is.na(minDate) & yearApproval >= minDate, TRUE, FALSE),
                postDate = ifelse(!is.na(minDate) & yearApproval < minDate, TRUE, FALSE)) %>%
    group_by(yearApproval) %>%
    summarise(c_humanTargets = n(),
              c_hasAnyGE = sum(hasAnyGE),
              c_priorDate = sum(priorDate), 
              c_postDate = sum(postDate)) %>%
    left_join(data_all_count) %>%
    mutate(c_nonhumanTargets = n - c_humanTargets,
            f_nonhumanTargets = c_nonhumanTargets/n * 100,
            c_noGE = c_humanTargets - c_hasAnyGE,
            f_noGE = c_noGE/n * 100,
            f_priorDate = c_priorDate/n * 100,
            f_postDate = c_postDate/n * 100)

data_humT_fraction <- data_humT_metadata %>%            
    select(yearApproval, f_noGE, f_nonhumanTargets, f_postDate, f_priorDate) %>%
    rename(without_GE = f_noGE,
           noH = f_nonhumanTargets,
           GE_noH_after = f_postDate,
           GE_noH_before = f_priorDate) %>%
    gather("type", "fraction", without_GE, noH, GE_noH_after, GE_noH_before) 

data_humT_count <- data_humT_metadata %>%            
    select(yearApproval, c_noGE, c_nonhumanTargets, c_postDate, c_priorDate) %>%
    rename(without_GE = c_noGE,
           noH = c_nonhumanTargets,
           GE_noH_after = c_postDate,
           GE_noH_before = c_priorDate) %>%
    gather("type", "count", without_GE, noH, GE_noH_after, GE_noH_before) 


data_humT_fraction %>% left_join(data_humT_count, by = c("yearApproval", "type")) %>%
    group_by(yearApproval) %>%
    mutate(label_ypos = 100 - cumsum(fraction) + fraction / 2) %>%
    ungroup() %>%
    mutate(yearApproval = as.factor(yearApproval), 
           type = as.factor(type)) 

custom_order <- c("noH", "without_GE", "GE_noH_after", "GE_noH_before")

data_to_plot <- data_humT_fraction %>% 
    mutate(type = fct_relevel(type, custom_order)) %>%
    arrange(type) %>%
    left_join(data_humT_count, by = c("yearApproval", "type")) %>%
    group_by(yearApproval) %>%
    mutate(label_ypos = 100 - (cumsum(fraction) - fraction / 2)) %>%
    ungroup() %>%
    mutate(yearApproval = as.factor(yearApproval), 
           type = as.factor(type))


# my_colors <- factor(c("#F39B7F99", "#8491B499", "#3C548899"))
my_colors <- c("#4a83af", "#6caad9","#f2b7ae", "#dddddd")


GE <- data_to_plot %>%
ggplot(aes(x=yearApproval, y=fraction, fill=fct_relevel(type, custom_order))) +
  geom_bar(stat="identity", width = 0.7) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0, 0), labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = my_colors,
                    breaks = c("GE_noH_before", "GE_noH_after", "without_GE", "noH"),
                    labels = c("Prior evidence", "No prior evidence", "No evidence", "Non-human target"),
                    name = "Genetic evidence: ") +
  geom_text(aes(x = yearApproval, y = label_ypos, label = ifelse(count>1, count, "")), size = 3.5, family=font, color = "#303030") +
  theme_minimal() +
  ylab("FDA approved drugs") +
  theme(axis.text.x = element_text(vjust = -0.6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.25, colour = "#737373"),
        axis.ticks = element_line(size = 0.25, color="#737373"),
        axis.ticks.length = unit(0.1, "cm"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        legend.position = "bottom",
        legend.text = element_text(family=font, size = 10),
        legend.title  = element_blank(),
        # = element_text(size = 10.5, face = "bold"),
        # legend.direction = "vertical",
        text = element_text(family=font, size = 12, color = "#4D4D4D"))

ggsave("./results/GE-Year_plot.png", 
        width = 6.5,
        height = 4.4,
        dpi = 600)