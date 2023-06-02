library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(epitools)
library(ggpubr)



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

data_humT_1 <- data_humT_metadata %>%            
    select(Year, without_GE_fraction, noH_fraction, GE_fraction_noH) %>%
    gather("type", "value", without_GE_fraction, noH_fraction, GE_fraction_noH) %>%
    mutate(Year = as.factor(Year), type = as.factor(type))

my_colors <- c("#F39B7F99", "#8491B499", "#3C548899")
# rev(RColorBrewer::brewer.pal(4, "Blues")[2:4])

ggplot(data = data_humT_1, aes(x=Year, y=value, fill=fct_reorder(type, value))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_colors, labels = c("Without evidence", "Non-human target", "With evidence")) +
  ylab("Approvals, %") +
  labs(fill='Evidence presence')

ggsave("./output/nature_plot/2013-2022_approvals_v2.pdf", 
        width = 6,
        height = 3)







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

colnames(table_data_P) <- c("Other", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("Other", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

colnames(table_data_O) <- c("Other", "Orphan")
rownames(table_data_O) <- c("no_GE", "has_GE")

colnames(table_data_S) <- c("Other", "Standard")
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
  odds_ratio_df <- data.frame(table_data = paste(rev(colnames(table_data)), collapse = "_vs_"),
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

ggplot(df, aes(x = estimate, y = table_data)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.2, 10)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x = "Odds Ratio (log scale)", y = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  geom_text(aes(label = paste0("p = ", ifelse(fisher_p_value < 0.001, format(signif(fisher_p_value, 2), scientific = TRUE), round(fisher_p_value, 3))), x = upper * 1.4),
            hjust = 0, size = 6)


ggsave("./data/nature_plot/OR_v1.png", 
        width = 10,
        height = 5)
