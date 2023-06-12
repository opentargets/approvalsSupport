library("tidyverse")

data1 <- read_csv("./data/2013-2022/2016-2022_approvals_v1.csv")
data2 <- read_csv("./data/2013-2022/2014-2015_approvals_v2.2.csv")
data3 <- read_csv("./data/2013-2022/2013_approvals_v2.csv")

local_approvals <- bind_rows(data1, data2, data3)

write.table(local_approvals, sep = ",", file = "./data/2013-2022/2013-2022_approvals_v1.csv", row.names = FALSE)


# Merge sells
df <- read_csv("./data/2013-2022/v3/2013-2022_approvals_GE_v3_in.csv")

df <- df %>%
  unite(Review_type_2, Review_type, Fast_track, sep = ", ", na.rm = TRUE)

write.table(df, sep = ",", file = "./data/2013-2022/v3/2013-2022_approvals_GE_v3.2_in.csv", row.names = FALSE)
