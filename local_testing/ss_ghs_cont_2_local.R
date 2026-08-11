## SS GHS CONTINUOUS -- LOCAL TEST, SINGLE SEED, MEDMIX CVs GIVEN

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)

setwd("C:/Users/01459189/OneDrive/phd/mixed-stratification") ## Update with WD
source("core_functions/function_calc_variance.R")

# ------------------------ TEST COMBINATION -------------------
max_clusters <- 15   # <-- upper bound only; KmeansSolution2() picks the best
seed <- 1            # <-- single seed for local testing

dir.create("OUTPUT", showWarnings = FALSE)

#------------------------ Data (seed-independent, loaded once) ---------------
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

# ------------------------ Single-seed run ----------
  
  ############## LOAD IN PARAMETERS FROM MEDMIX_GA OUTPUT (this seed only)
  filename <- paste0(
    "OUTPUT/medmix_ga_cont_", "seed", seed, "_results_LOCALTEST.Rdata"
  )
  load(filename)  # brings in `store_out` from the medmix run for this seed
  
  ssize <- sum(store_out$n)
  cv_y <- unlist(store_out$cv_y)
  cv_b <- unlist(store_out$cv_b)
  cv_g <- unlist(store_out$cv_g)
  cv <- cbind(t(cv_y), cv_b, cv_g)
  
  df_cv <- as.data.frame(cv)
  colnames(df_cv) <- c("cv_y_1", "cv_y_2", "cv_b_1", "cv_g_1")
  
  ## Format data nicely:
  vars <- character()
  cfs <- character()
  sf <- ghs
  cv_list <- list(df_cv)
  error <- as.matrix(cv_list[[1]], nrow = 1)
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
  frame3 <- buildFrameDF(
    df = sf,
    id = "id",
    X = c("X1", "X2", "X3", "X4"),
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
  
  #------------------------ KmeansSolution2 initial solution --------------------
  set.seed(seed)
  kmean <- KmeansSolution2(
    frame = frame3,
    errors = error,
    maxclusters = max_clusters
  )
  
  nstrat <- tapply(kmean$suggestions, kmean$domainvalue, FUN = function(x) length(unique(x)))
  
  ## Prepare suggestion for optimisation
  sugg <- prepareSuggestion(kmean, frame3, nstrat)
  
  #------------------------ Run continuous GA (single seed, single run) -------
  start.time <- Sys.time()
  proc.time.start <- proc.time()
  
  solution <- optimStrata(
    method = "continuous",
    errors = error,
    framesamp = frame3,
    nStrata = nstrat,
    suggestions = sugg,
    iter = 2000,
    pops = 50,
    elitism_rate = 0.2,
    mut_chance = 0.7
  )
  
  proc_time <- proc.time() - proc.time.start
  end.time <- Sys.time()
  time_taken <- difftime(end.time, start.time, units = "mins")
  
  strataStructure <- summaryStrata(
    solution$framenew,
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
  
  # ----------------------- Store Results (single run, this seed) --------------
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
    cv_y_fpc = fitness$cv_k_fpc,
    cv_b_fpc = fitness$cv_b_fpc,
    cv_g_fpc = fitness$cv_g_fpc,
    proc_time = proc_time,
    used_initial_solution = TRUE
  )
  
  outfile <- paste0(
    "OUTPUT/ss_ghs_cont2_seed", seed, "_results_LOCALTEST.Rdata"
  )
 # save(store, file = outfile)
  
  list(seed = seed, status = "ok", outfile = outfile,
       n_strata_realized = n_strata_realized, time_mins = as.numeric(time_taken))
  
}