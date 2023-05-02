library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
# install.packages("epitools")
library(epitools)

data_in <- read_csv("./output/2013-2022_approvals_v3.csv")
data  <- subset(data_in, targetIds != "NA" & has_GE != "NA" & TA != "Oncology")

# data$P_O <- grepl("P", data$Review_type) & grepl("O", data$Review_type)
# table_data_P_O  <- table(data$has_GE, data$P_O)
# odds_ratio_P_O <- odds_ratio_df(table_data_P_O)


# data$P_noO <- grepl("P(?!.*O)", data$Review_type, perl = TRUE)
# table_data_P_noO  <- table(data$has_GE, data$P_noO)
# odds_ratio_P_noO <- odds_ratio_df(table_data_P_noO)


# Create boolean columns with Review status and modalities
data$Priority <- grepl("P", data$Review_type)
data$Breakthrough <- grepl("B", data$Review_type)
# data$Accelerated <- grepl("A", data$Review_type)
data$Orphan <- grepl("O", data$Review_type)
data$SmallMol <- grepl("Small molecule", data$drugType)
data$Ab <- grepl("Antibody", data$drugType)
data$Protein <- grepl("Protein|Enzyme", data$drugType)
data$Oligonucl <- grepl("Oligonucleotide", data$drugType)
data$Standard <- grepl("S", data$Review_type)

# write.table(data, sep = ",", file = "./output/2013-2022_approvals_v4.csv", row.names = FALSE)


# create a 2x2 contingency table with "has_GE" values as rows and the review status values as columns.
table_data_P <- table(data$has_GE, data$Priority)
table_data_B <- table(data$has_GE, data$Breakthrough)
# table_data_A <- table(data$has_GE, data$Accelerated) # Odds ration is not applicable
table_data_O <- table(data$has_GE, data$Orphan)
table_data_S <- table(data$has_GE, data$Standard)
table_data_SmallMol <- table(data$has_GE, data$SmallMol)
table_data_Ab <- table(data$has_GE, data$Ab)
table_data_Protein <- table(data$has_GE, data$Protein)
table_data_Oligonucl <- table(data$has_GE, data$Oligonucl)

# table_data_SmallMol
# table_data_Ab
# table_data_Protein
# table_data_Oligonucl

colnames(table_data_P) <- c("Other", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("Other", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

# colnames(table_data_A) <- c("Other", "Accelerated")
# rownames(table_data_A) <- c("no_GE", "has_GE")

colnames(table_data_O) <- c("Other", "Orphan")
rownames(table_data_O) <- c("no_GE", "has_GE")

colnames(table_data_S) <- c("Other", "Standard")
rownames(table_data_S) <- c("no_GE", "has_GE")

colnames(table_data_SmallMol) <- c("Other", "SmallMol")
rownames(table_data_SmallMol) <- c("no_GE", "has_GE")

colnames(table_data_Ab) <- c("Other", "Ab")
rownames(table_data_Ab) <- c("no_GE", "has_GE")

colnames(table_data_Protein) <- c("Other", "Protein")
rownames(table_data_Protein) <- c("no_GE", "has_GE")

colnames(table_data_Oligonucl) <- c("Other", "Oligonucl")
rownames(table_data_Oligonucl) <- c("no_GE", "has_GE")


# dimnames(table_data_P) <- list(has_GE = c(TRUE, FALSE), Priority = c(TRUE, FALSE))
# dimnames(table_data_B) <- list(has_GE = c(TRUE, FALSE), Breakthrough = c(TRUE, FALSE))
# dimnames(table_data_A) <- list(has_GE = c(TRUE, FALSE), Accelerated = c(TRUE, FALSE))
# dimnames(table_data_O) <- list(has_GE = c(TRUE, FALSE), Orphan = c(TRUE, FALSE))

# fisher_result <- fisher.test(table_data_A)

# # Create a vector of correction factors
# correction_factor <- ifelse(table_data_A < 5, 1, 0)

# # Apply Yates' correction to the contingency table_data_A
# table_data_A <- table_data_A + correction_factor


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
# odds_ratio_4 <- odds_ratio_df(table_data_A)
odds_ratio_5 <- odds_ratio_df(table_data_SmallMol)
odds_ratio_6 <- odds_ratio_df(table_data_Ab)
odds_ratio_7 <- odds_ratio_df(table_data_Protein)
odds_ratio_8 <- odds_ratio_df(table_data_Oligonucl)
odds_ratio_9 <- odds_ratio_df(table_data_S)


# Combine results into one data frame
all_odds_ratios <- rbind(odds_ratio_1, odds_ratio_2, odds_ratio_3, odds_ratio_9)

# View results
all_odds_ratios

write.table(all_odds_ratios, sep = ",", file = "./output/2013-2022_approvals_v4.4_nonco.csv", row.names = FALSE)







all_odds_ratios <- read_csv("./output/2013-2022_approvals_v4.4_nonco.csv")

df <- all_odds_ratios[order(all_odds_ratios$estimate), ]



# Create a factor variable for table_data to use in the y-axis
df$table_data <- factor(df$table_data, levels = df$table_data)

ggplot(df, aes(x = estimate, y = table_data)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(limits = c(0.2, 10)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x = "Odds Ratio (log scale)", y = "", title = "2013-2022 approvals (non oncology only)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  geom_text(aes(label = paste0("p = ", round(fisher_p_value, 3)), x = upper + 0.3),
            hjust = 0, size = 6)


ggsave("./output/OR_v4.4.png", 
        width = 10,
        height = 4)
