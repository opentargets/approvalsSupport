library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")

data <- read_csv("./data/nature_plot/2013-2022_approvals_v2.csv") %>%
        distinct(Drug_name, .keep_all = TRUE) 

# selected_cols <- data %>%
#    select(Drug_name, noTarget)

local_approvals <- read_csv("./data/nature_plot/2013-2022_approvals_v3.csv")

approvals_nHt <- merge(x = local_approvals, y = data[ , c("Drug_name", "noTarget")], by = "Drug_name", all.x=TRUE)


write.table(approvals_nHt, sep = ",", file = "./data/nature_plot/2013-2022_approvals_v3.2_nHt.csv", row.names = FALSE)



