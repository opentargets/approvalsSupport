library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
# install.packages("epitools")
library(epitools)

data_in <- read_csv("./output/2013-2022_approvals_v2_hasGE.csv")
data  <- subset(data_in, targetIds != "NA" & has_GE != "NA")


# Create boolean columns with Review status
data$Priority <- grepl("P", data$Review_type)
data$Breakthrough <- grepl("B", data$Review_type)
data$Accelerated <- grepl("A", data$Review_type)
data$Orphan <- grepl("O", data$Review_type)

# write.table(data, sep = ",", file = "./output/2013-2022_approvals_v2_review.csv", row.names = FALSE)


# create a 2x2 contingency table with "has_GE" values as rows and the review status values as columns.
table_data_P <- table(data$has_GE, data$Priority)
table_data_B <- table(data$has_GE, data$Breakthrough)
table_data_A <- table(data$has_GE, data$Accelerated)
table_data_O <- table(data$has_GE, data$Orphan)

colnames(table_data_P) <- c("Other", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("Other", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

colnames(table_data_A) <- c("Other", "Accelerated")
rownames(table_data_A) <- c("no_GE", "has_GE")

colnames(table_data_O) <- c("Other", "Orphan")
rownames(table_data_O) <- c("no_GE", "has_GE")


# dimnames(table_data_P) <- list(has_GE = c(TRUE, FALSE), Priority = c(TRUE, FALSE))
# dimnames(table_data_B) <- list(has_GE = c(TRUE, FALSE), Breakthrough = c(TRUE, FALSE))
# dimnames(table_data_A) <- list(has_GE = c(TRUE, FALSE), Accelerated = c(TRUE, FALSE))
# dimnames(table_data_O) <- list(has_GE = c(TRUE, FALSE), Orphan = c(TRUE, FALSE))

# fisher_result <- fisher.test(table_data)

# # Create a vector of correction factors
correction_factor <- ifelse(table_data_A < 5, 1, 0)

# Apply Yates' correction to the contingency table_data_A
table_corrected_A <- table_data_A + correction_factor


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
  odds_ratio_df <- data.frame(table_data = paste(colnames(table_data), collapse = "_vs_"),
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
odds_ratio_4 <- odds_ratio_df(table_corrected_A)

# Combine results into one data frame
all_odds_ratios <- rbind(odds_ratio_1, odds_ratio_2, odds_ratio_3, odds_ratio_4)

# # View results
# all_odds_ratios


all_odds_ratios <- all_odds_ratios[order(all_odds_ratios$estimate), ]

# Create a new column with the y-axis position of each group
all_odds_ratios$y <- seq(nrow(all_odds_ratios))

# Plot the forest plot
ggplot(all_odds_ratios, aes(x = estimate, y = table_data)) +
  geom_point(size = 4, color = "blue") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "blue") +
  xlim(c(0, 16)) +
  scale_y_reverse(breaks = seq(nrow(df)), labels = all_odds_ratios$group) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey", size = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")


ggsave("./output/OR_v1.png", 
        width = 6,
        height = 3)
