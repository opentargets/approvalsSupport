library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(epitools)
library(ggpubr)

font = "Helvetica"

data_in <- read_csv("./results/2013-2022_approvals_preGE.csv")

data  <- subset(data_in, hasAnyGE != "NA")  

# Create boolean columns with Review status type and disease type
data$Accelerated <- grepl("A", data$reviewType)
data$Priority <- grepl("P", data$reviewType)
data$Breakthrough <- grepl("B", data$reviewType)
data$Fast_track <- grepl("F", data$reviewType)
data$Expedited <- !(grepl("^S$|^S,O$", data$reviewType))
data$Serious <- grepl("Serious", data$diseaseClass)

# Create 2x2 contingency tables with "hasAnyGE" values as rows and the review status or disease type as columns
table_data_A <- table(data$hasAnyGE, data$Accelerated)
table_data_P <- table(data$hasAnyGE, data$Priority)
table_data_B <- table(data$hasAnyGE, data$Breakthrough)
table_data_F <- table(data$hasAnyGE, data$Fast_track)
table_data_E <- table(data$hasAnyGE, data$Expedited)
table_data_Serious <- table(data$hasAnyGE, data$Serious)

# Rename 2x2 contingency tables
colnames(table_data_A) <- c("", "Accelerated")
rownames(table_data_A) <- c("no_GE", "has_GE")

colnames(table_data_P) <- c("", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

colnames(table_data_F) <- c("", "Fast track")
rownames(table_data_F) <- c("no_GE", "has_GE")

colnames(table_data_E) <- c("", "Expedited")
rownames(table_data_E) <- c("no_GE", "has_GE")

colnames(table_data_Serious) <- c("", "Serious conditions")
rownames(table_data_Serious) <- c("no_GE", "has_GE")

# # Merge 2x2 contingency tables
# list_of_matrix_to_merge <- lapply(list(
#   table_data_A, 
#   table_data_P,
#   table_data_B,
#   table_data_F,
#   table_data_E,
#   table_data_Serious), as.data.frame.matrix)
# merged_tibble <- do.call(cbind, list_of_matrix_to_merge)

# write.csv(merged_tibble, file = "./results/OR_2x2.csv", row.names = TRUE)


# Calculate odds ratios
odds_ratio_df <- function(table_data) {
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
odds_ratio_A <- odds_ratio_df(table_data_A)
odds_ratio_P <- odds_ratio_df(table_data_P)
odds_ratio_B <- odds_ratio_df(table_data_B)
odds_ratio_F <- odds_ratio_df(table_data_F)
odds_ratio_E <- odds_ratio_df(table_data_E)
odds_ratio_SD <- odds_ratio_df(table_data_Serious)

# Combine results into one data frame
all_odds_ratios <- rbind(odds_ratio_A,
                          odds_ratio_P,
                          odds_ratio_B,
                          odds_ratio_F,
                          odds_ratio_E,
                          odds_ratio_SD)
df <- all_odds_ratios[order(all_odds_ratios$estimate), ]

# write.table(df, sep = ",", file = "./results/OR_values.csv", row.names = FALSE)

# Create a factor variable for table_data to use in the y-axis
df$table_data <- factor(df$table_data, levels = df$table_data)

# Plot ORs
OR <- ggplot(df, aes(x = estimate, y = table_data)) +
  geom_point(colour = "#4D4D4D") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.4, height = 0.1, colour = "#4D4D4D") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#4d4d4d95") +
  scale_x_log10(limits = c(0.7, 60), labels = function(x) ifelse(x > 1, round(x), x)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x = "Odds Ratio (log scale)", y = "") +
  theme(axis.text.x = element_text(vjust = -0.6),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
        axis.title.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#ffffff", colour = NA),
        axis.line = element_line(size = 0.25, colour = "#737373"),
        axis.ticks = element_line(size = 0.25, color="#737373"),
        axis.ticks.length = unit(0.1, "cm"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        legend.position = "bottom",
        # legend.direction = "vertical",
        text=element_text(family=font, size = 12, color = "#4D4D4D"),
        plot.margin = margin(t = 4, r = 0, b = 0, l = 0)) +
geom_text(aes(label = paste0("p = ", ifelse(fisher_p_value < 0.001, format(signif(fisher_p_value, 2), scientific = TRUE), round(fisher_p_value, 3))),
              x = upper * 1.3), hjust = 0, size = 3.1, color = "#4D4D4D")
              

ggsave("./results/OR_plot.png", 
        width = 3,
        height = 3)
