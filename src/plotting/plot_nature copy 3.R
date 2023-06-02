library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(epitools)
library(ggpubr)

font = "Helvetica"

## PLOT 1 (GE decade)
data <- read_csv("./data/nature_plot/2013-2022_approvals_v2.csv")

# Calculation of overall approval numbers per year
data_all_count <- data %>%
    select(Drug_name_original, Year) %>%
    distinct() %>%
    group_by(Year) %>%
    summarise(n = n()) %>%
    ungroup()

# Creating a table with % of has_GE, has_no_GE and non-human targets
data_humT_metadata <- data %>% filter(noTarget == FALSE) %>%
    mutate(score_bool = score > 0, 
           interactionAssociation_bool = replace_na(interactionAssociation, FALSE),
           phenotypeAssociation_bool = replace_na(phenotypeAssociation, FALSE),
          #  has_GE = interactionAssociation_bool|score_bool) %>%
           has_GE = interactionAssociation_bool|phenotypeAssociation_bool|score_bool) %>%
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

data_humT_fraction <- data_humT_metadata %>%            
    select(Year, without_GE_fraction, noH_fraction, GE_fraction_noH) %>%
    rename(without_GE = without_GE_fraction,
           noH = noH_fraction,
           GE_noH = GE_fraction_noH) %>%
    gather("type", "fraction", without_GE, noH, GE_noH) 

data_humT_count <- data_humT_metadata %>%            
    select(Year, without_GE_count, noH, count_GE) %>%
    rename(without_GE = without_GE_count,
           GE_noH = count_GE) %>%
    gather("type", "count", without_GE, noH, GE_noH) 


data_humT_fraction %>% left_join(data_humT_count, by = c("Year", "type")) %>%
    group_by(Year) %>%
    mutate(label_ypos = 100 - cumsum(fraction) + fraction / 2) %>%
    ungroup() %>%
    mutate(Year = as.factor(Year), 
           type = as.factor(type)) %>%
    filter(Year == 2013)

data_to_plot <- data_humT_fraction %>% left_join(data_humT_count, by = c("Year", "type")) %>%
    group_by(Year) %>%
    mutate(label_ypos = 100 - (cumsum(fraction) - fraction / 2)) %>%
    ungroup() %>%
    mutate(Year = as.factor(Year), 
           type = as.factor(type))

# my_colors <- factor(c("#F39B7F99", "#8491B499", "#3C548899"))
my_colors <- c("#6caad9", "#b8ccda", "#f2b7ae")
# color_mapping <- c("With evidence" = "6caad9", "Without evidence" = "f2b7ae", "Non-human target" = "b8ccda")


GE <- data_to_plot %>%
ggplot(aes(x=Year, y=fraction, fill=fct_reorder(type, fraction))) +
  geom_bar(stat="identity", width = 0.7) +
  scale_y_continuous(limits = c(0,100.01), expand = c(0, 0), labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = my_colors,
                    breaks = c("GE_noH", "noH", "without_GE"),
                    labels = c("With evidence", "Non-human target", "Without evidence"),
                    name = "Evidence presence: ") +
  geom_text(aes(x = Year, y = label_ypos, label = count), size = 3.5, family=font, color = "#303030") +
  theme_minimal() +
  ylab("Approvals") +
  theme(axis.text.x = element_text(vjust = -0.6),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.25, colour = "#737373"),
        axis.ticks = element_line(size = 0.25, color="#737373"),
        axis.ticks.length = unit(0.1, "cm"),
        plot.background = element_rect(fill = "#F7F7F7", color = NA),
        legend.position = "bottom",
        legend.text = element_text(family=font, size = 10),
        legend.title = element_text(size = 10.5, face = "bold"),
        # legend.direction = "vertical",
        text = element_text(family=font, size = 12, color = "#4D4D4D"))

# ggsave("./output/nature_plot/2013-2022_approvals_v2.pdf", 
#         width = 5,
#         height = 5)


### PLOT 2 (OR)
data_in <- read_csv("./data/nature_plot/2013-2022_approvals_v3.csv")
data  <- subset(data_in, targetIds != "NA" & has_GE != "NA")  %>% 
        distinct(Brand_name, .keep_all = TRUE) 



# Create boolean columns with Review status and modalities
data$Priority <- grepl("P", data$Review_type)
data$Breakthrough <- grepl("B", data$Review_type)
data$Orphan <- grepl("O", data$Review_type)
data$Standard <- grepl("S", data$Review_type)
# data$Accelerated <- grepl("A", data$Review_type)

# create a 2x2 contingency table with "has_GE" values as rows and the review status values as columns.
table_data_P <- table(data$has_GE, data$Priority)
table_data_B <- table(data$has_GE, data$Breakthrough)
table_data_O <- table(data$has_GE, data$Orphan)
table_data_S <- table(data$has_GE, data$Standard)

# table_data_A <- table(data$has_GE, data$Accelerated) # Odds ration is not applicable

colnames(table_data_P) <- c("", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

colnames(table_data_O) <- c("", "Orphan")
rownames(table_data_O) <- c("no_GE", "has_GE")

colnames(table_data_S) <- c("", "Standard")
rownames(table_data_S) <- c("no_GE", "has_GE")

odds_ratio_df <- function(table_data) {
  # Calculate odds ratio
  odds_ratio_result <- oddsratio(table_data)
  
  # Extract estimate, lower, and upper values
  odds_ratio_estimate <- odds_ratio_result$measure[2, 1]
  odds_ratio_lower <- odds_ratio_result$measure[2, 2]
  odds_ratio_upper <- odds_ratio_result$measure[2, 3]
  
  # Extract midp.exact, fisher.exact, and chi.square p-values
  midp_p_value <- odds_ratio_result$p.value[2, 1]
  fisher_p_value <- odds_ratio_result$p.value[2, 2]
  chi_square_p_value <- odds_ratio_result$p.value[2, 3]
  
  # Create data frame with odds ratio results
  odds_ratio_df <- data.frame(table_data = paste(rev(colnames(table_data)), collapse = ""),
                               estimate = odds_ratio_estimate,
                               lower = odds_ratio_lower,
                               upper = odds_ratio_upper,
                               midp_p_value = midp_p_value,
                               fisher_p_value = fisher_p_value,
                               chi_square_p_value = chi_square_p_value,
                               stringsAsFactors = FALSE)
  
  # Return data frame
  return(odds_ratio_df)
}

# Calculate odds ratio and extract values for each table
odds_ratio_1 <- odds_ratio_df(table_data_P)
odds_ratio_2 <- odds_ratio_df(table_data_B)
odds_ratio_3 <- odds_ratio_df(table_data_O)
odds_ratio_4 <- odds_ratio_df(table_data_S)

# Combine results into one data frame
all_odds_ratios <- rbind(odds_ratio_1,
                        odds_ratio_2,
                        odds_ratio_3,
                        odds_ratio_4)

df <- all_odds_ratios[order(all_odds_ratios$estimate), ]

# Create a factor variable for table_data to use in the y-axis
df$table_data <- factor(df$table_data, levels = df$table_data)

OR <- ggplot(df, aes(x = estimate, y = table_data)) +
  geom_point(colour = "#4D4D4D") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.4, height = 0.1, colour = "#4D4D4D") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#4d4d4d95") +
  scale_x_log10(limits = c(0.2, 14), labels = function(x) ifelse(x > 1, round(x), x)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x = "Odds Ratio (log scale)", y = "") +
#   theme_bw() +
  theme(axis.text.x = element_text(vjust = -0.6),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#F7F7F7", colour = NA),
        axis.line = element_line(size = 0.25, colour = "#737373"),
        axis.ticks = element_line(size = 0.25, color="#737373"),
        axis.ticks.length = unit(0.1, "cm"),
        plot.background = element_rect(fill = "#F7F7F7", color = NA),
        legend.position = "bottom",
        # legend.direction = "vertical",
        text=element_text(family=font, size = 12, color = "#4D4D4D"),
        plot.margin = margin(t = 0, r = 10, b = 0, l = -40)) +
geom_text(aes(label = paste0("p = ", ifelse(fisher_p_value < 0.001, format(signif(fisher_p_value, 2), scientific = TRUE), round(fisher_p_value, 3))),
              x = upper * 1.3), hjust = 0, size = 3.1, color = "#4D4D4D") 

# ggsave("./output/nature_plot/OR_v1.png", 
#         width = 10,
#         height = 5)



library(patchwork)


combined <- GE + OR + 
        plot_layout(guides = "collect", widths = c(1.8, 1)) & 
        theme(legend.position = "bottom") 

ggsave("./output/nature_plot/GE_OR_v5.pdf", 
        width = 8,
        height = 3.9)
