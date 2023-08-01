library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(epitools)
library(ggpubr)

font = "Helvetica"

### PLOT 2 (OR)
data_in <- read_csv("./output/v3/2013-2022_approvals_GE_v3.4_out.csv")

data  <- subset(data_in, has_GE != "NA")  

# Create boolean columns with Review status and modalities
data$Accelerated <- grepl("A", data$Review_type)
data$Priority <- grepl("P", data$Review_type)
data$Breakthrough <- grepl("B", data$Review_type)
data$Orphan <- grepl("O", data$Review_type)
data$Standard <- grepl("S", data$Review_type) #below standard is everything non-expedited
data$nExpedited <- grepl("^S$|^S, O$", data$Review_type)
data$Fast_track <- grepl("F", data$Review_type)
data$Expedited <- !(grepl("^S$|^S, O$", data$Review_type))
data$Oncology <- grepl("Oncology", data$TA)
data$nOncology <- grepl("Other", data$TA)


# create a 2x2 contingency table with "has_GE" values as rows and the review status values as columns.
table_data_A <- table(data$has_GE, data$Accelerated)
table_data_P <- table(data$has_GE, data$Priority)
table_data_B <- table(data$has_GE, data$Breakthrough)
table_data_O <- table(data$has_GE, data$Orphan)
table_data_S <- table(data$has_GE, data$Standard)
table_data_F <- table(data$has_GE, data$Fast_track)
table_data_E <- table(data$has_GE, data$Expedited)
table_data_nE <- table(data$has_GE, data$nExpedited)
table_data_onco <- table(data$has_GE, data$Oncology)
table_data_nonco <- table(data$has_GE, data$nOncology)


colnames(table_data_A) <- c("", "Accelerated")
rownames(table_data_A) <- c("no_GE", "has_GE")

colnames(table_data_P) <- c("", "Priority")
rownames(table_data_P) <- c("no_GE", "has_GE")

colnames(table_data_B) <- c("", "Breakthrough")
rownames(table_data_B) <- c("no_GE", "has_GE")

colnames(table_data_O) <- c("", "Orphan")
rownames(table_data_O) <- c("no_GE", "has_GE")

colnames(table_data_S) <- c("", "Standard")
rownames(table_data_S) <- c("no_GE", "has_GE")

colnames(table_data_F) <- c("", "Fast track")
rownames(table_data_F) <- c("no_GE", "has_GE")

colnames(table_data_E) <- c("", "Expedited")
rownames(table_data_E) <- c("no_GE", "has_GE")

colnames(table_data_nE) <- c("", "non-Expedited")
rownames(table_data_nE) <- c("no_GE", "has_GE")

colnames(table_data_onco) <- c("", "Oncology")
rownames(table_data_onco) <- c("no_GE", "has_GE")

colnames(table_data_nonco) <- c("", "non-Oncology")
rownames(table_data_nonco) <- c("no_GE", "has_GE")

# lapply for each of cont. table matrix function, merge them by using do.call (to call cbind iteratively)
list_of_matrix_to_merge <- lapply(list(
  table_data_A, 
  table_data_P,
  table_data_B,
  table_data_O,
  table_data_S,
  table_data_F,
  table_data_E,
  table_data_nE,
  table_data_onco,
  table_data_nonco), as.data.frame.matrix)
merged_tibble <- do.call(cbind, list_of_matrix_to_merge)

write.csv(merged_tibble, file = "./output/v3/OR/OR_v3.4_all_2x2.csv", row.names = TRUE)

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
odds_ratio_0 <- odds_ratio_df(table_data_A)
odds_ratio_1 <- odds_ratio_df(table_data_P)
odds_ratio_2 <- odds_ratio_df(table_data_B)
odds_ratio_3 <- odds_ratio_df(table_data_O)
odds_ratio_4 <- odds_ratio_df(table_data_S)
odds_ratio_5 <- odds_ratio_df(table_data_F)
odds_ratio_6 <- odds_ratio_df(table_data_E)
odds_ratio_7 <- odds_ratio_df(table_data_nE)
odds_ratio_8 <- odds_ratio_df(table_data_onco)
odds_ratio_9 <- odds_ratio_df(table_data_nonco)

# Combine results into one data frame
all_odds_ratios <- rbind(
                        odds_ratio_0,
                        odds_ratio_1,
                        odds_ratio_2,
                        odds_ratio_3,
                        odds_ratio_4,
                        odds_ratio_5,
                        odds_ratio_6,
                        odds_ratio_7,
                        odds_ratio_8,
                        odds_ratio_9)

df <- all_odds_ratios[order(all_odds_ratios$estimate), ]

write.table(df, sep = ",", file = "./output/v3/OR/OR_v3.4_all.csv", row.names = FALSE)

# Create a factor variable for table_data to use in the y-axis
df$table_data <- factor(df$table_data, levels = df$table_data)

OR <- ggplot(df, aes(x = estimate, y = table_data)) +
  geom_point(colour = "#4D4D4D") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.4, height = 0.1, colour = "#4D4D4D") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#4d4d4d95") +
  scale_x_log10(limits = c(0.05, 130), labels = function(x) ifelse(x > 1, round(x), x)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  labs(x = "Odds Ratio (log scale)", y = "") +
#   theme_bw() +
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
              

ggsave("./output/v3/OR/OR_v3.4_all_tm.png", 
        width = 7,
        height = 4)
