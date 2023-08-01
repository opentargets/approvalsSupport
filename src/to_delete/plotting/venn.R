library("tidyverse")
library("ggvenn")
library(ggplot2)

data_in <- read_csv("./output/v3/2013-2022_approvals_GE_v3.4_out.csv")
df  <- subset(data_in, has_GE != "NA")

df$onco <- grepl("Oncology", df$TA)
df$Expedited <- !(grepl("^S$|^S, O$", df$Review_type))

# write.table(df, sep = ",", file = "./output/v3/2013-2022_approvals_GE_v3.4_out.csv", row.names = FALSE)


# P_count <- sum(df$Priority)
# B_count <- sum(df$Breakthrough)
# A_count <- sum(df$Accelerated)
# F_count <- sum(df$Fast_Track)
# O_count <- sum(df$Orphan)
# GE_count <- sum(df$has_GE)
onco_count <- sum(df$onco)
expedited_count <- sum(df$Expedited)

# # Expedited + oncology + GE
# expedited_GE <- sum(df$Expedited & df$has_GE)
# onco_GE <- sum(df$Expedited & df$has_GE)
expedited_onco <- sum(df$Expedited & df$onco)

# GE <- as.vector(df$has_GE)
# expedited <- as.vector(df$Expedited)
# onco <- as.vector(df$onco)


# Expedited + oncology + GE
data_venn_1 <- df[, c("onco", "Expedited", "has_GE")]
ggplot(data_venn_1,
       aes(A = onco, B = Expedited, C = has_GE)) + 
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
data_venn_3 <- df[, c("Expedited", "Orphan", "has_GE")]
ggplot(data_venn_3,
       aes(A = Expedited, B = Orphan, C = has_GE)) + 
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
data_venn_1_2 <- df[, c("onco", "Expedited")]
ggplot(data_venn_1_2,
       aes(A = onco, B = Expedited)) + 
  theme_void() +
  geom_venn()

ggsave("./output/v3/Exp_onc_venn_1.png", 
        width = 5,
        height = 5)


# Fast_Track + Orphan 
data_venn_4 <- df[, c("Fast_Track", "Orphan")]
ggplot(data_venn_4,
       aes(A = Fast_Track, B = Orphan)) + 
  theme_void() +
  geom_venn()

ggsave("./output/F_O_v1.png", 
        width = 3,
        height = 3)
