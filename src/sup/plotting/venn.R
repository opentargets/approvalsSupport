library("tidyverse")
library("ggvenn")
library(ggplot2)

data_in <- read_csv("./data/nature_plot/2013-2022_approvals_v3.2_ft.csv")
df  <- subset(data_in, targetIds != "NA" & has_GE != "NA")  %>% 
        distinct(Brand_name, .keep_all = TRUE) 

df$onco <- grepl("Oncology", df$TA)
df$all_noS <- !(grepl("^S$|^S, O$", df$Review_type) & !df$Fast_Track)

write.table(df, sep = ",", file = "./output/2013-2022_approvals_v3_3.csv", row.names = FALSE)


# P_count <- sum(df$Priority)
# B_count <- sum(df$Breakthrough)
# A_count <- sum(df$Accelerated)
# F_count <- sum(df$Fast_Track)
# O_count <- sum(df$Orphan)
# GE_count <- sum(df$has_GE)
onco_count <- sum(df$onco)
expedited_count <- sum(df$all_noS)

# # Expedited + oncology + GE
# expedited_GE <- sum(df$all_noS & df$has_GE)
# onco_GE <- sum(df$all_noS & df$has_GE)
expedited_onco <- sum(df$all_noS & df$onco)

# GE <- as.vector(df$has_GE)
# expedited <- as.vector(df$all_noS)
# onco <- as.vector(df$onco)


# Expedited + oncology + GE
data_venn_1 <- df[, c("onco", "all_noS", "has_GE")]
ggplot(data_venn_1,
       aes(A = onco, B = all_noS, C = has_GE)) + 
  theme_void() +
  geom_venn()

ggsave("./output/Exp_onc_GE_v1.png", 
        width = 5.5,
        height = 5.5)

# Orphan + oncology + GE
data_venn_2 <- df[, c("onco", "Orphan", "has_GE")]
ggplot(data_venn_2,
       aes(A = onco, B = Orphan, C = has_GE)) + 
  theme_void() +
  geom_venn()

ggsave("./output/O_onc_GE_v1.png", 
        width = 5.5,
        height = 5.5)


# Orphan + Expedited + GE
data_venn_3 <- df[, c("all_noS", "Orphan", "has_GE")]
ggplot(data_venn_3,
       aes(A = all_noS, B = Orphan, C = has_GE)) + 
  theme_void() +
  geom_venn()

ggsave("./output/Exp_O_GE_v1.png", 
        width = 5.5,
        height = 5.5)



# Expedited types
data_venn <- df[, c("Priority", "Breakthrough", "Accelerated", "Fast_Track")]
ggplot(data_venn,
       aes(A = Priority, B = Breakthrough, C = Accelerated, D = Fast_Track)) +
  geom_venn() +
  theme_void()

ggsave("./output/Exp_types_v1.png", 
        width = 10,
        height = 8)



# Expedited + oncology 
data_venn_1_2 <- df[, c("onco", "all_noS")]
ggplot(data_venn_1_2,
       aes(A = onco, B = all_noS)) + 
  theme_void() +
  geom_venn()

ggsave("./output/Exp_onc_v1.png", 
        width = 3,
        height = 3)


# Fast_Track + Orphan 
data_venn_4 <- df[, c("Fast_Track", "Orphan")]
ggplot(data_venn_4,
       aes(A = Fast_Track, B = Orphan)) + 
  theme_void() +
  geom_venn()

ggsave("./output/F_O_v1.png", 
        width = 3,
        height = 3)
