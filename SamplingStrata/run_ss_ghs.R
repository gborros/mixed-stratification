## Run File - Parallel version across num_strata x seed grid
## Adds tryCatch-wrapped kmeans initial solution for the continuous method

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)
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

ncores <- min(40, nrow(param_grid)) # don't request more workers than tasks

# ------------------------ Cluster setup -----------------------------------
# PSOCK cluster: each worker is a fresh R process, so packages and sourced
# files must be loaded on each worker via clusterEvalQ.

cl <- makeCluster(ncores)
clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(cluster)
  library(SamplingStrata)
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
  
  type1 <- list(numeric = 1, numeric = 2, symm = 3, factor = 4)
  type_name1 <- c("numeric", "numeric", "symm", "factor")
  
  datasets <- c("ghs")
  dfs <- list(ghs)
  type <- list(type1)
  type_name <- list(type_name1)
  
  ## CVs
  load("CVs/cv_input_app.RData")
  
  df_cv <- df_cv_min %>%
    filter(
      type == "original",
      strata == num_strata,
      N == ssize
    ) %>%
    select(cv_y_1, cv_y_2, cv_b_1, cv_g_1)
  
  vars <- character()
  cfs <- character()
  sf <- ghs
  cv <- list(df_cv)
  error <- as.matrix(cv[[1]], nrow = 1)
  data_type <- type[[1]]
  data_name <- type_name[[1]]
  
  for (z in 1:ncol(sf)) {
    var <- paste0("X", z)
    vars <- c(vars, var)
    cf <- paste0("CV", z)
    cfs[z] <- paste0(cf)
  }
  
  colnames(sf) <- vars
  colnames(error) <- cfs
  
  sf$id <- 1:nrow(sf) ## adding identifier
  sf$dom <- 1 ## adding domain
  
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
  
  #------------------------ Kmeans initial solution (optional) ------------
  set.seed(seed)
  initial_solution <- tryCatch(
    {
      sol <- KmeansSolution2(
        frame = frame3,
        errors = error,
        nstrat = num_strata
      )
      
      n_groups <- length(unique(unlist(sol)))
      
      if (n_groups != num_strata) {
        message(sprintf(
          "[num_strata=%d, seed=%d] KmeansSolution2 produced %d distinct strata (expected %d) — discarding suggestion.",
          num_strata, seed, n_groups, num_strata
        ))
        NULL
      } else {
        sol
      }
    },
    error = function(e) {
      message(sprintf(
        "[num_strata=%d, seed=%d] KmeansSolution2 failed: %s — proceeding without initial solution.",
        num_strata, seed, conditionMessage(e)
      ))
      NULL
    }
  )
  
  #------------------------ Run GGA (single seed per task) -----------------
  
  set.seed(seed)
  start.time <- Sys.time()
  proc.time.start <- proc.time()
  
  optimStrata_args <- list(
    method = "continuous",
    errors = error,
    framesamp = frame3,
    iter = 2000,
    pops = 50,
    nStrata = num_strata,
    elitism_rate = 0.2,
    mut_chance = 0.7
  )
  
  if (!is.null(initial_solution)) {
    optimStrata_args$suggestions <- initial_solution
  }
  
  used_initial_flag <- !is.null(initial_solution)
  
  solution <- tryCatch(
    do.call(optimStrata, optimStrata_args),
    error = function(e) {
      if (!is.null(initial_solution) &&
          grepl("suggestions not compatible with nStrata", conditionMessage(e))) {
        message(sprintf(
          "[num_strata=%d, seed=%d] optimStrata rejected suggestions (%s) — retrying without them.",
          num_strata, seed, conditionMessage(e)
        ))
        used_initial_flag <<- FALSE
        optimStrata_args$suggestions <- NULL
        do.call(optimStrata, optimStrata_args)
      } else {
        stop(e)
      }
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
  n_strata_realized <- length(n) # actual number of strata this run produced
  
  df_sol <- solution$framenew
  df_sol$strata <- df_sol$STRATO
  
  fitness <- calculate_variance(
    df = df_sol, df_sol$strata,
    vars = vars, n = n, type = data_type, type_name = data_name, seed = seed
  )
  
  # ----------------------- Store Results (single run per task) --------------
  # Each task corresponds to exactly one num_strata x seed combination,
  
  store <- list(
    num_strata = num_strata,
    seed = seed,
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
    n_strata_realized = n_strata_realized,
    used_initial_solution = used_initial_flag
  )
  
  filename <- paste0(
    "OUTPUT/ss_", datasets, "_", num_strata, "strata_seed", seed, "_", ssize, "n", "_results.Rdata"
  )
  
  save(store, file = filename)
  
  filename # return the filename so the foreach result list shows what was written
}

stopCluster(cl)