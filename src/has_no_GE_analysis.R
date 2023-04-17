library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
# install.packages("epitools")
library(epitools)

data_in <- read_csv("./output/2013-2022_approvals_v1.2_mod.csv")
data  <- subset(data_in, targetIds != "NA" & has_GE != "NA")

# write.table(data, sep = ",", file = "./output/test1.csv", row.names = FALSE)

# Standard approvals

data$has_S <- grepl("S", data$Review_type)

# create a 2x2 contingency table with "has_GE" values as columns and the "has_S" values as rows.
table_data <- table(data$has_S, data$has_GE)

# calculate the odds ratio using the oddsratio function
odds_ratio <- oddsratio(table_data)

# print the odds ratio (The output will be the odds ratio along with its confidence interval, 
# which tells you how much more likely it is for an individual to have the outcome 
# (has_S = TRUE) if they have the exposure (has_GE = TRUE) compared to if they do not have the exposure.)
odds_ratio



# Orphan approvals
data$has_O <- grepl("O", data$Review_type)

# create a 2x2 contingency table with "has_GE" values as rows and the "has_S" values as columns.
table_data <- table(data$has_O, data$has_GE)

# calculate the odds ratio using the oddsratio function
odds_ratio <- oddsratio(table_data)

# print the odds ratio (The output will be the odds ratio along with its confidence interval, 
# which tells you how much more likely it is for an individual to have the outcome 
# (has_S = TRUE) if they have the exposure (has_GE = TRUE) compared to if they do not have the exposure.)
odds_ratio



# Breakthrough approvals
data$has_B <- grepl("B", data$Review_type)

# create a 2x2 contingency table with "has_GE" values as rows and the "has_S" values as columns.
table_data <- table(data$has_B, data$has_GE)

# calculate the odds ratio using the oddsratio function
odds_ratio <- oddsratio(table_data)

# print the odds ratio (The output will be the odds ratio along with its confidence interval, 
# which tells you how much more likely it is for an individual to have the outcome 
# (has_S = TRUE) if they have the exposure (has_GE = TRUE) compared to if they do not have the exposure.)
odds_ratio




data_in <- read_csv("./output/2013-2022_approvals_v1.2_mod.csv")
data  <- subset(data_in, targetIds != "NA" & has_GE != "NA")

# Small molecules
data$drugType_s <- grepl("Small molecule", data$drugType)


# create a 2x2 contingency table with "has_GE" values as columns and the "has_S" values as rows.
table_data <- table(data$drugType_s, data$has_GE)

# calculate the odds ratio using the oddsratio function
odds_ratio <- oddsratio(table_data)

# print the odds ratio (The output will be the odds ratio along with its confidence interval, 
# which tells you how much more likely it is for an individual to have the outcome 
# (has_S = TRUE) if they have the exposure (has_GE = TRUE) compared to if they do not have the exposure.)
odds_ratio
