## SS GHS ATOMIC LOCAL
## v2: final number of strata is whatever KmeansSolution() suggests, not forced.
## max_clusters below is only an UPPER BOUND for KmeansSolution's internal search -
## it is not the final strata count. Everything downstream (prints, filename) now
## refers to n_strata_realized, the actual outcome, not this cap.

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)

setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

# ------------------------ SET YOUR TEST COMBINATION HERE -------------------
max_clusters <- 15   # <-- upper bound only; KmeansSolution() picks the best <= this. Assumption: 15, adjust as you like.
seed         <- 3    # <-- change this

# ------------------------ Source local functions ---------------------------
source("core_functions/medmix_gasampsi.R")
source("core_functions/function_calc_variance.R")

#------------------------ Data ---------------------------------
cat(">> Loading GHS data...\n")
load("multivar_datasets/ghs_2024_mixed.RData")
ghs <- as.data.frame(df)
cat(sprintf("   Raw dims: %d rows x %d cols\n", nrow(ghs), ncol(ghs)))

ghs <- ghs[ghs$fin_reqinc != 9999999, ]
cat(sprintf("   After dropping fin_reqinc == 9999999: %d rows\n", nrow(ghs)))

ghs$fin_reqinc <- as.numeric(ghs$fin_reqinc)
ghs <- ghs[, -4]
ghs <- ghs[, -3]
ghs <- ghs[, -3]
cat(sprintf("   After column drops: %d cols -> %s\n", ncol(ghs), paste(colnames(ghs), collapse = ", ")))
ghs$fin_reqinc_cat <- var.bin(ghs$fin_reqinc, 15)

ghs$hwl_status <- as.factor(ghs$hwl_status)
ghs$geotype <- as.factor(ghs$geotype)
ghs$head_age <- as.integer(ghs$head_age)
ghs$head_age_cat <- var.bin(ghs$head_age, 15)

ghs$hwl_status <- as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
ghs$geotype <- as.integer(ghs$geotype == 1) ## urban

cat("   Post-recode summary:\n")
print(summary(ghs))

type1 <- list(numeric = 1, numeric = 2, symm = 3, factor = 4)
type_name1 <- c("numeric", "numeric", "symm", "factor")
nvars <- length(type_name1)

########## LOAD INPUTS FROM MEDMIX
filename <- paste0(
  "OUTPUT/medmix_ga_atomic_", "seed", seed, "_results_LOCALTEST.Rdata"
)
load(filename)  # brings in `store` from the medmix run for this seed

ssize <- sum(store_out$n)
cv_y <- unlist(store_out$cv_y)
cv_b <- unlist(store_out$cv_b)
cv_g <- unlist(store_out$cv_g)
cv <- cbind(t(cv_y), cv_b, cv_g)
strata <- store_out$strata

df_cv <- as.data.frame(cv)
colnames(df_cv) <- c("cv_y_1", "cv_y_2", "cv_b_1", "cv_g_1")

vars <- character()
cfs <- character()
sf <- ghs
cv <- list(df_cv)
error <- as.matrix(cv[[1]], nrow = 1)
data_type <- type1
data_name <- type_name1

for (z in 1:nvars) {
  var <- paste0("X", z)
  vars <- c(vars, var)
  cf <- paste0("CV", z)
  cfs[z] <- paste0(cf)
}

colnames(sf) <- c(vars, "X1_cat", "X2_cat")
colnames(error) <- cfs

sf$id <- 1:nrow(sf) ## adding identifier
sf$dom <- 1 ## adding domain

cat("\n>> Building sampling frame with buildFrameDF()...\n")
frame3 <- buildFrameDF(
  df = sf,
  id = "id",
  X = c("X1_cat", "X2_cat", "X3", "X4"),
  Y = c("X1", "X2", "X3", "X4"),
  domainvalue = "dom"
)
cat(sprintf("   Frame built: %d rows x %d cols\n", nrow(frame3), ncol(frame3)))
str(frame3)

## specify error thresholds:
ndom <- length(unique(sf$dom))
error <- as.data.frame(list(
  DOM = rep("DOM1", ndom),
  error,
  domainvalue = c(1)
))
cat("\n>> Error/CV threshold table passed to optimStrata:\n")
print(error)

## ADD MY OWN STRATA

frame3$STRATO <- strata  # length must equal nrow(frame3), one label per unit

strata_df <- frame3 %>%
  group_by(domainvalue, STRATO) %>%
  summarise(
    N   = n(),
    M1  = mean(Y1), S1 = sd(Y1),
    M2  = mean(Y2), S2 = sd(Y2),
    M3  = mean(Y3), S3 = sd(Y3),
    M4  = mean(Y4), S4 = sd(Y4),
    .groups = "drop"
  ) %>%
  mutate(
    COST = 1,
    CENS = 0,
    DOM1 = domainvalue
  ) %>%
  dplyr::select(-domainvalue)

strata_df <- as.data.frame(strata_df)
str(strata_df)

alloc <- bethel(
  stratif = strata_df,
  errors  = error,          # the same error/cv threshold data frame you already built
  minnumstrat = 2,
  printa = TRUE
)

n = alloc