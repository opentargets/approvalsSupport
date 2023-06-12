library("cowplot")
library("ggsci")
library(forcats)
library(ggplot2)
library("tidyverse")
library(ggpubr)

data <- read_csv("./data/nature_plot/2013-2022_approvals_v3.2_ft.csv")

data_all_count <- data %>%
    select(Drug_name_original, Year) %>%
    distinct() %>%
    group_by(Year) %>%
    summarise(n = n()) %>%
    ungroup()

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

data_humT_2 <- data_humT_metadata %>%            
    select(Year, without_GE_count, noH, count_GE) %>%
    gather("type", "value", without_GE_count, noH, count_GE) %>%
    mutate(Year = as.factor(Year), 
    type = as.factor(type),
    label_ypos=cumsum(value))


my_colors <- c("#f3aca2", "#dadada", "#84b6dd")
# rev(RColorBrewer::brewer.pal(4, "Blues")[2:4])

ggplot(data = data_humT_1, aes(x=Year, y=value, fill=fct_reorder(type, value))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_colors, labels = c("Without evidence", "Non-human target", "With evidence")) +
  ylab("Approvals, %") +
  labs(fill='Evidence presence')

ggsave("./output/GE_2013-2022_approvals_v1.2_1.pdf", 
        width = 6,
        height = 3)


ggplot(data = data_humT_2, aes(x=Year, y=value, fill=fct_reorder(type, value))) +
  geom_bar(stat="identity") +
  theme_minimal() +
  scale_fill_manual(values = my_colors, labels = c("Without evidence", "Non-human target", "With evidence")) +
  ylab("Approvals numbers") +
  labs(fill='Evidence presence')


ggsave("./output/GE_2013-2022_approvals_v1.2_2.pdf", 
        width = 6,
        height = 3)



# GE linear regression (GE_fraction_noH)
subset_df <- subset(data_humT_1, type == "GE_fraction_noH")
# subset_df <- subset(subset_df, Year != 2013)


fit = lm(value ~ Year, data = subset_df)

# subset_df$Year <- factor(subset_df$Year, levels = sort(levels(subset_df$Year)))
# numeric_years <- as.numeric(subset_df$Year)

# numeric_years <- as.numeric(as.character(Year))

ggplot(subset_df ,aes(x = as.integer(Year), y = as.numeric(value))) +
  xlab("Year") + 
  ylab("GE, %") +
  geom_smooth(method='lm', color = "#84b6dd") +
  theme_minimal(base_size = 14) +
  geom_point() +
  geom_text(label = as.integer(subset_df$value), nudge_y = 2.5) +
  scale_x_continuous(limits = c(1, 10)) +
  # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
  #                    "Intercept =",signif(fit$coef[[1]],5 ),
  #                    " Slope =",signif(fit$coef[[2]], 5),
  #                    " P =",signif(summary(fit)$coef[2,4], 5))) 
  stat_cor(label.x = 2, label.y = 78) +
  stat_regline_equation(label.x = 2, label.y = 75)

ggsave("./output/GE_2013-2022_regr_v1.pdf", 
        width = 6,
        height = 3)




# GE linear regression (without_GE_fraction, noH_fraction, GE_fraction_noH)
subset_df <- subset(data_humT_1, type == "without_GE_fraction")
# subset_df <- subset(subset_df, Year != 2013)

fit = lm(value ~ Year, data = subset_df)

# subset_df$Year <- factor(subset_df$Year, levels = sort(levels(subset_df$Year)))
# numeric_years <- as.numeric(subset_df$Year)

# numeric_years <- as.numeric(as.character(Year))

ggplot(subset_df ,aes(x = as.integer(Year), y = as.numeric(value))) +
  xlab("Year") + 
  ylab("GE, %") +
  geom_smooth(method='lm', color = "#f3aca2") +
  theme_minimal(base_size = 14) +
  geom_point() +
  geom_text(label = as.integer(subset_df$value), nudge_y = 2.5) +
  scale_x_continuous(limits = c(1, 10)) +
  # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
  #                    "Intercept =",signif(fit$coef[[1]],5 ),
  #                    " Slope =",signif(fit$coef[[2]], 5),
  #                    " P =",signif(summary(fit)$coef[2,4], 5))) 
  stat_cor(label.x = 4, label.y = 39) +
  stat_regline_equation(label.x = 4, label.y = 42)

ggsave("./output/NoGE_2013-2022_regr_v1.pdf", 
        width = 6,
        height = 3)



# GE linear regression (noH_fraction)
subset_df <- subset(data_humT_1, type == "noH_fraction")
# subset_df <- subset(subset_df, Year != 2013)

fit = lm(value ~ Year, data = subset_df)

# subset_df$Year <- factor(subset_df$Year, levels = sort(levels(subset_df$Year)))
# numeric_years <- as.numeric(subset_df$Year)

# numeric_years <- as.numeric(as.character(Year))

ggplot(subset_df ,aes(x = as.integer(Year), y = as.numeric(value))) +
  xlab("Year") + 
  ylab("GE, %") +
  geom_smooth(method='lm', color = "#4d4d4d") +
  theme_minimal(base_size = 14) +
  geom_point() +
  geom_text(label = as.integer(subset_df$value), nudge_y = 2.5) +
  scale_x_continuous(limits = c(1, 10)) +
  # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
  #                    "Intercept =",signif(fit$coef[[1]],5 ),
  #                    " Slope =",signif(fit$coef[[2]], 5),
  #                    " P =",signif(summary(fit)$coef[2,4], 5))) 
  stat_cor(label.x = 4, label.y = 39) +
  stat_regline_equation(label.x = 4, label.y = 42)

ggsave("./output/NoHT_2013-2022_regr_v1.pdf", 
        width = 6,
        height = 3)
