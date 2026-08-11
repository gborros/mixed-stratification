# Local Testing of MED-GA method using GGA-BC (atomic) inputs

library(dplyr)
library(haven)
library(cluster)
library(doParallel)
library(foreach)

setwd("C:/Users/01459189/OneDrive/phd/mixed-stratification") ## Update with WD
source("core_functions/function_calc_variance.R")
source("core_functions/medmix_gasampsi.R")
  #------------------------ Data ---------------------------------

  ## GHS
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

  type1 <- list(numeric = 1, numeric = 2, symm = 3, factor = 4)
  type_name1 <- c("numeric", "numeric", "symm", "factor")

  datasets <- c("ghs")
  dfs <- list(ghs)
  type <- list(type1)
  type_name <- list(type_name1)
  
  ############## LOAD IN PARAMETERS FROM SAMPLINGSTRATA OUTPUT
  
  seed <- 1 ## Update with seed wanting to test
  
  filename <- paste0(
    "OUTPUT/ss_ghs_atomic1_", "seed", seed, "_results_LOCALTEST.Rdata"
  )
  
  load(filename)
  
  ssize = sum(store$n)
  num_strata = store$n_strata_realized
  
  ########## CONTINUE WITH DATA SET-UP

  vars <- character()
  sf <- dfs[[1]]
  data_type <- type[[1]]
  data_name <- type_name[[1]]

  for (z in 1:ncol(sf)) {
    var <- paste0("X", z)
    vars <- c(vars, var)
  }

  colnames(sf) <- vars
  df_in <- sf

  #------------------------ Run GGA (single seed per task) -----------------

  start.time <- Sys.time()
  proc.time.start <- proc.time()

  result <- run_medmix_gasampsi(
    seed = seed, df = df_in, pop_size = 50, num_strata = num_strata,
    num_generations = 2000, mutation_rate = 0.7,
    vars = vars, ssize = ssize, elitism_rate = 0.2,
    crossover_rate = 0.7, type = data_type, type_name = data_name
  )

  proc_time <- proc.time() - proc.time.start
  end.time <- Sys.time()
  time_taken <- difftime(end.time, start.time, units = "mins")

  strata <- result$best_solution

  # ----------------------- Store Results (single run per task) --------------
  # Each task corresponds to exactly one num_strata x seed combination,

  store <- list(
    num_strata = num_strata,
    seed = seed,
    n = result$best_n,
    strata = strata,
    time = time_taken,
    fitness = result$best_fitness,
    fitness_log = result$fitness_log,
    deff_y = result$best_deff_y,
    deff_b = result$best_deff_b,
    deff_g = result$best_deff_g,
    cv_y = result$best_cv_y,
    cv_b = result$best_cv_b,
    cv_g = result$best_cv_g,
    proc_time = proc_time
  )

  filename <- paste0(
    "OUTPUT/medmix_ga_atomic_", "seed", seed, "_results_LOCALTEST.Rdata"
  )

 # save(store, file = filename)
