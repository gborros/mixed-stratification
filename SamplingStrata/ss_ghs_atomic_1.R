## SS GHS ATOMIC -- 10 SEEDS IN PARALLEL, STANDARD CVS GIVEN

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)
library(doParallel)
library(foreach)

setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")
source("core_functions/function_calc_variance.R")

# ------------------------ TEST COMBINATION -------------------
max_clusters <- 15   # <-- upper bound only; KmeansSolution() picks the best
seeds        <- 1:10 # <-- 10 seeds run in parallel

#------------------------ Data (seed-independent, built ONCE) ---------------
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

## BIN continuous vars:
set.seed(1234)
ghs$fin_reqinc_cat <- var.bin(ghs$fin_reqinc, 15)
ghs$head_age_cat <- var.bin(ghs$head_age, 15)

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

colnames(sf) <- c(vars, "X1_cat", "X2_cat")
colnames(error) <- cfs

sf$id <- 1:nrow(sf) ## adding identifier
sf$dom <- 1 ## adding domain

## Make sample frame for SamplingStrata (seed-independent -- built once):
frame3 <- buildFrameDF(
  df = sf,
  id = "id",
  X = c("X1_cat", "X2_cat", "X3", "X4"),
  Y = c("X1", "X2", "X3", "X4"),
  domainvalue = "dom"
)

## specify error thresholds:
ndom <- length(unique(sf$dom))
error <- as.data.frame(list(
  DOM = rep("DOM1", ndom),
  error,
  domainvalue = c(1)
))

#------------------------ Build atomic strata (buildStrataDF, seed-independent) ---------------
strata_atomic <- buildStrataDF(frame3, progress = FALSE)

#------------------------ Register parallel backend --------------------------
n_cores <- 10
cl <- makeCluster(n_cores, outfile = "OUTPUT/parallel_log.txt")
registerDoParallel(cl)

dir.create("OUTPUT", showWarnings = FALSE)

#------------------------ Run atomic GA across seeds in parallel -------------
results <- foreach(
  seed = seeds,
  .packages = c("SamplingStrata", "dplyr"),
  .errorhandling = "pass"
) %dopar% {
  
  ## each worker needs its own copy of helper functions/env
  source("core_functions/function_calc_variance.R")
  
  set.seed(seed)
  kmean <- KmeansSolution(
    strata = strata_atomic,
    errors = error,
    maxclusters = max_clusters
  )
  
  nstrat <- tapply(kmean$suggestions, kmean$domainvalue, FUN = function(x) length(unique(x)))
  
  start.time <- Sys.time()
  proc.time.start <- proc.time()
  
  solution <- tryCatch(
    optimStrata(
      method = "atomic",
      errors = error,
      framesamp = frame3,
      nStrata = nstrat,
      suggestions = kmean,
      iter = 2000,
      pops = 50,
      elitism_rate = 0.2,
      mut_chance = 0.7
    ),
    error = function(e) {
      message(sprintf(
        "[seed=%d] optimStrata(method='atomic') failed: %s",
        seed, conditionMessage(e)
      ))
      stop(e)
    }
  )
  
  proc_time <- proc.time() - proc.time.start
  end.time <- Sys.time()
  time_taken <- difftime(end.time, start.time, units = "mins")
  
  strataStructure <- summaryStrata(solution$framenew,
                                   solution$aggr_strata,
                                   progress = FALSE
  )
  
  n <- as.vector(strataStructure$Allocation)
  n_strata_realized <- length(n)
  
  df_sol <- solution$framenew
  df_sol$strata <- df_sol$LABEL
  
  fitness <- calculate_variance(
    df = df_sol, df_sol$strata,
    vars = vars, n = n, type = data_type, type_name = data_name, seed = seed
  )
  
  # ----------------------- Store Results (per seed) --------------
  store <- list(
    max_clusters = max_clusters,
    seed = seed,
    method = "atomic",
    n_atomic_strata = nrow(strata_atomic),
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
    cv_y_fpc = fitness$cv_k_fpc,
    cv_b_fpc = fitness$cv_b_fpc,
    cv_g_fpc = fitness$cv_g_fpc,
    proc_time = proc_time,
    used_initial_solution = TRUE
  )
  
  filename <- paste0(
    "OUTPUT/ss_ghs_atomic1_", "seed", seed, "_results_LOCALTEST.Rdata"
  )
  save(store, file = filename)
  
  store
}

stopCluster(cl)
