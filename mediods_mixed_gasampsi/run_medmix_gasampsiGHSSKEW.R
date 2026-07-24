## Run File - Parallel version across num_strata x seed grid

library(dplyr)
library(haven)
library(cluster)
library(doParallel)
library(foreach)

#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

# ------------------------ Grid setup -------------------------------------
# Each row of param_grid = one (num_strata, seed) combination = one parallel task.
# 4 strata values x 10 seeds = 40 tasks, matching 40 cores.
# Change seed_values if you want a different number of tasks.

num_strata_values <- c(3, 4, 5, 6)
seed_values <- 1:10

param_grid <- expand.grid(
  num_strata = num_strata_values,
  seed = seed_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

ncores <- min(12, nrow(param_grid)) # don't request more workers than tasks

# ------------------------ Cluster setup -----------------------------------
# PSOCK cluster: each worker is a fresh R process, so packages and sourced
# files must be loaded on each worker via clusterEvalQ.

cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(cluster)
  source("core_functions/medmix_gasampsi.R")
  source("core_functions/function_calc_variance.R")
})
registerDoParallel(cl)

# ------------------------ Parallel loop over param_grid --------------------

foreach(row_i = seq_len(nrow(param_grid)), .errorhandling = "pass") %dopar% {

  num_strata <- param_grid$num_strata[row_i]
  seed <- param_grid$seed[row_i]

  ssize <- 2000

  #------------------------ Data ---------------------------------

  ## GHS
  load("multivar_datasets/ghs_2024_mixed.RData")
  ghs <- as.data.frame(df)
  ghs <- ghs[ghs$fin_reqinc!=9999999,]
  ghs$fin_reqinc <- as.numeric(ghs$fin_reqinc)
  ghs <- ghs[,-4]
  ghs <- ghs[,-3]
  ghs <- ghs[,-3]
  ghs <- ghs[,-1]
  
  ghs$hwl_status <- as.factor(ghs$hwl_status)
  ghs$geotype <- as.factor(ghs$geotype)
  
  ghs$hwl_status <-as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
  ghs$geotype <-as.integer(ghs$geotype==1) ## urban
  
  
  type1 <- list(numeric=1, symm=2, factor=3)
  type_name1 <- c("numeric", "symm", "factor")
  
  datasets <- c("ghs_skew")
  dfs <- list(ghs)
  type <- list(type1)
  type_name <- list(type_name1)

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
    "OUTPUT/medmix_ga_", datasets, "_", num_strata, "strata_seed", seed, "_", ssize, "n", "_results.Rdata"
  )

  save(store, file = filename)

  filename # return the filename so the foreach result list shows what was written
}

stopCluster(cl)
