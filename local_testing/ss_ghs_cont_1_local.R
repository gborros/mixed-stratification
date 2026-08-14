## SS GHS CONTINUOUS LOCAL -- FIRST ROUND WITH STANDARD 0.05 CVS GIVEN

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)

setwd("C:/Users/01459189/OneDrive/phd/mixed-stratification") ## Update with WD
source("core_functions/function_calc_variance.R")

# ------------------------ TEST COMBINATION -------------------
max_clusters <- 15   # <-- upper bound only; KmeansSolution2() picks the best
seed         <- 1    # <-- change this

#------------------------ Data ---------------------------------
load("multivar_datasets/ghs_2024_mixed.RData")
ghs <- as.data.frame(df)
ghs <- ghs[ghs$fin_reqinc != 9999999, ]
ghs$fin_reqinc <- as.numeric(ghs$fin_reqinc)
ghs <- ghs[, -4]
ghs <- ghs[, -3]
ghs <- ghs[, -3]
ghs$hwl_status <- as.factor(ghs$hwl_status)
ghs$geotype <- as.factor(ghs$geotype)
ghs$head_age <- as.integer(ghs$head_age)
ghs$hwl_status <- as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
ghs$geotype <- as.integer(ghs$geotype == 1) ## urban

ghs <- ghs[,-5]

print(summary(ghs))

## Inputs for local functions:
type1 <- list(numeric = 1, numeric = 2, symm = 3, factor = 4)
type_name1 <- c("numeric", "numeric", "symm", "factor")
nvars <- length(type_name1)

## INPUT OWN CVS
cv <- c(0.05, 0.05, 0.05, 0.05)
df_cv <- as.data.frame(t(cv))
colnames(df_cv) <- c("cv_y_1", "cv_y_2", "cv_b_1", "cv_g_1")

## Format data nicely:
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

colnames(sf) <- c(vars)
colnames(error) <- cfs

sf$id <- 1:nrow(sf) ## adding identifier
sf$dom <- 1 ## adding domain


## Make sample frame for SamplingStrata:
cat("\n>> Building sampling frame with buildFrameDF()...\n")
frame3 <- buildFrameDF(
  df = sf,
  id = "id",
  X = c("X1", "X2", "X3", "X4"),
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

#------------------------ Build atomic strata (buildStrataDF) ---------------
# cat("\n>> Building atomic strata with buildStrataDF()...\n")
# strata_atomic <- buildStrataDF(frame3, progress = TRUE)
# cat(sprintf("\n   buildStrataDF() returned %d atomic strata (frame has %d units).\n",
#             nrow(strata_atomic), nrow(frame3)))

#------------------------ KmeansSolution initial solution --------------------
cat(sprintf("\n>> Running KmeansSolution() on atomic strata (maxclusters = %d, upper bound only)...\n", max_clusters))
set.seed(seed)
kmean <- KmeansSolution2(
  frame = frame3,
  errors = error,
  maxclusters = max_clusters
)
cat("   KmeansSolution() output (head):\n")
print(head(kmean))

nstrat <- tapply(kmean$suggestions, kmean$domainvalue, FUN = function(x) length(unique(x)))
cat("\n>> Final nStrata per domain, as suggested by KmeansSolution (this is what gets used):\n")
print(nstrat)

# Prepare suggestion for optimization
sugg <- prepareSuggestion(kmean,frame3,nstrat)

#------------------------ Run atomic GA (single seed, single run) -----------
cat("\n>> Running optimStrata(method = 'atomic')...\n")
set.seed(seed)
start.time <- Sys.time()
proc.time.start <- proc.time()

solution <- tryCatch(
  optimStrata(
    method = "continuous",
    errors = error,
    framesamp = frame3,
    nStrata = nstrat,
    suggestions = sugg,
    iter = 2000,
    pops = 50,
    elitism_rate = 0.2,
    mut_chance = 0.7
  ),
  error = function(e) {
    message(sprintf(
      "[seed=%d] optimStrata(method='continuous') failed: %s",
      seed, conditionMessage(e)
    ))
    stop(e)
  }
)

proc_time <- proc.time() - proc.time.start
end.time <- Sys.time()
time_taken <- difftime(end.time, start.time, units = "mins")
cat(sprintf("\n>> optimStrata (continuous) finished in %.2f minutes.\n", as.numeric(time_taken)))

cat("\n>> Structure of solution object:\n")
str(solution, max.level = 1)

strataStructure <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress = FALSE
)
cat("\n>> summaryStrata() output:\n")
print(strataStructure)

n <- as.vector(strataStructure$Allocation)
n_strata_realized <- length(n) # actual number of strata this run produced
cat(sprintf("\n>> KmeansSolution suggested nStrata (per domain) = %s | Realized strata after GA = %d\n",
            paste(nstrat, collapse = ","), n_strata_realized))
cat(">> Allocation vector (n):\n")
print(n)
cat(sprintf(">> Total allocated sample size: %d\n", sum(n)))

df_sol <- solution$framenew
df_sol$strata <- df_sol$LABEL

cat("\n>> Strata frequency table (from df_sol$strata):\n")
print(table(df_sol$strata))

cat("\n>> Calculating variance / CV diagnostics via calculate_variance()...\n")
fitness <- calculate_variance(
  df = df_sol, df_sol$strata,
  vars = vars, n = n, type = data_type, type_name = data_name, seed = seed
)

cat("\n>> Fitness / CV results:\n")
cat(sprintf("   deff (overall) : %s\n", fitness$deff))
cat(sprintf("   deff_y         : %s\n", paste(fitness$deff_y, collapse = ", ")))
cat(sprintf("   deff_b         : %s\n", paste(fitness$deff_b, collapse = ", ")))
cat(sprintf("   deff_g         : %s\n", paste(fitness$deff_g, collapse = ", ")))
cat(sprintf("   cv_y (cv_k)    : %s\n", paste(fitness$cv_k, collapse = ", ")))
cat(sprintf("   cv_b           : %s\n", paste(fitness$cv_b, collapse = ", ")))
cat(sprintf("   cv_g           : %s\n", paste(fitness$cv_g, collapse = ", ")))

# ----------------------- Store Results (single run) --------------
store <- list(
  max_clusters = max_clusters,
  seed = seed,
  method = "continuous",
  nstrat_suggested = nstrat,
  n_strata_realized = n_strata_realized,
  df_sol = df_sol,
  n = n,
  strata = df_sol$strata,
  time = time_taken,
  fitness = fitness$deff,
  deff_y = fitness$deff_y,
  deff_b = fitness$deff_b,
  deff_g = fitness$deff_g,
  cv_y = fitness$cv_k,
  cv_b = fitness$cv_b,
  cv_g = fitness$cv_g,
  proc_time = proc_time,
  used_initial_solution = TRUE
)

filename <- paste0(
  "OUTPUT/ss_ghs_cont1_", "seed", seed, "_results_LOCALTEST.Rdata"
)

#save(store, file=filename)