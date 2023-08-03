library(dplyr)
library(arsenal)

# df_old <- read_csv("results/2013-2022_approvals_GE_prec_cur.csv") #%>% select(brandDrugName, hasAnyGE)
# df <- read_csv("results/2013-2022_approvals_GE_final.csv") #%>% select(brandDrugName, hasAnyGE)


df_old <- read_csv("results/2013-2022_approvals_GE_by_source.csv") %>%
    distinct() %>%
    arrange(yearApproval, )
df <- read_csv("results/2013-2022_approvals_GE_map.csv") #%>% select(brandDrugName, hasAnyGE)


# differences <- anti_join(df_old, df, by = c("brandDrugName", "hasAnyGE"))

differences <- anti_join(df_old, df, by = c("originalDrugName", "yearApproval", "datasourceId", "drugName","interactionAssociation","phenotypeAssociation","score"))
# Print the result
print(differences, n =100)


# differences <- anti_join(df, df_old, by = c("brandDrugName", "hasAnyGE"))

# # Print the result
# print(differences)




