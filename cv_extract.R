## CV Extraction - Mixed

options(scipen = 999)
library(ggplot2)
library(cluster)
library(haven)
library(dplyr)
library(knitr)
library(ggthemes)

setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed/ANALYSIS26")

load("C:/Users/01459189/OneDrive/phd/hpc - mixed/WIP26/med/med_all_res.RData")

res$dataset <- sub("^d_", "", res$dataset)

num_map <- c(
  mix0 = 1,
  mix1 = 2,
  mix2 = 3,
  mix3 = 2,
  mix4 = 3,
  mix5 = 4,
  mix6 = 4
)

res$num_var <- num_map[res$dataset]

df_cv <- data.frame(dataset = res$dataset, strata = res$strata, N = res$N, deff = res$deff,
                    cv_y_1 = res$cv_y_1, cv_y_2 = res$cv_y_2, 
                    cv_b_1 = res$cv_b_1, cv_b_2 = res$cv_b_2, cv_g_1 = res$cv_g_1)


df_cv_min <- df_cv %>%
  group_by(dataset, strata, N) %>%
  slice_min(order_by = deff, n = 1, with_ties = FALSE) %>%
  ungroup()

save(df_cv_min, file = "C:/Users/01459189/OneDrive/phd/hpc - mixed/CVs/cv_input.RData")