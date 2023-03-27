library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)

options(dplyr.width = Inf)
data <- read_csv("./output/2018-2022_approvals_v1.csv")

data_all_count <- data %>%
    select(Drug_name_original, Year) %>%
    distinct() %>%
    group_by(Year) %>%
    summarise(n = n()) %>%
    ungroup()data_humT_metadata

data_humT_metadata <- data %>% filter(noTarget == FALSE) %>%
    mutate(score_bool = score > 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE),
           has_GE = interactionAssociation_bool|score_bool) %>%
        #    has_GE = interactionAssociation_bool|phenotypeAssociation_bool|score_bool) %>%
    group_by(datasourceId, Drug_name_original, Year) %>%
    summarise(has_GE = any(has_GE)) %>% 
    ungroup() %>%
    group_by(Drug_name_original, Year) %>%
    summarise(has_GE = any(has_GE)) %>% 
    ungroup() %>%
    group_by(Year) %>%
    summarise(n_filtered = n(), count_GE = sum(has_GE)) %>%
    left_join(data_all_count) %>%
    mutate(GE_fraction_noH = count_GE / n * 100,
            noH = n - n_filtered,
            without_GE_count = n_filtered - count_GE,
            noH_fraction = noH / n * 100,
            without_GE_fraction = without_GE_count / n * 100) 

data_humT_1 <- data_humT_metadata %>%            
    select(Year, without_GE_fraction, noH_fraction, GE_fraction_noH) %>%
    gather("type", "value", without_GE_fraction, noH_fraction, GE_fraction_noH) %>%
    mutate(Year = as.factor(Year), type = as.factor(type))

data_humT_2 <- data_humT_metadata %>%            
    select(Year, without_GE_count, noH, count_GE) %>%
    gather("type", "value", without_GE_count, noH, count_GE) %>%
    mutate(Year = as.factor(Year), 
    type = as.factor(type),
    label_ypos=cumsum(value))


my_colors <- c("#f3aca2", "#dadada", "#84b6dd")
# rev(RColorBrewer::brewer.pal(4, "Blues")[2:4])

ggplot(data = data_humT_1, aes(x=Year, y=value, fill=fct_reorder(type, value))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_colors, labels = c("Without evidence", "Non-human target", "With evidence")) +
  ylab("Approvals, %") +
  labs(fill='Evidence presence')

ggsave("./output/direct_evidence_interact_1.pdf", 
        width = 5,
        height = 3)


ggplot(data = data_humT_2, aes(x=Year, y=value, fill=fct_reorder(type, value))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_colors, labels = c("Without evidence", "Non-human target", "With evidence")) +
  ylab("Approvals numbers") +
  labs(fill='Evidence presence')


ggsave("./output/direct_evidence_interact_2.pdf", 
        width = 5,
        height = 3)


