## BETHEL GHS ATOMIC -- 10 SEEDS IN PARALLEL, USES MEDMIX-SUGGESTED STRATA

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)
library(doParallel)
library(foreach)

#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")
source("core_functions/medmix_gasampsi.R")
source("core_functions/function_calc_variance.R")

# ------------------------ TEST COMBINATION -------------------
seeds <- 1:10 # <-- 10 seeds run in parallel

#------------------------ Data (seed-independent, built ONCE) ---------------
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
ghs <- ghs[,-5]
cat(sprintf("   After column drops: %d cols -> %s\n", ncol(ghs), paste(colnames(ghs), collapse = ", ")))

set.seed(1234) ## set seed for consistent binning
ghs$fin_reqinc_cat <- var.bin(ghs$fin_reqinc, 15)
ghs$hwl_status <- as.factor(ghs$hwl_status)
ghs$geotype <- as.factor(ghs$geotype)
ghs$head_age <- as.integer(ghs$head_age)
ghs$head_age_cat <- var.bin(ghs$head_age, 15)
ghs$hwl_status <- as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
ghs$geotype <- as.integer(ghs$geotype == 1) ## urban

cat("   Post-recode summary:\n")
print(summary(ghs))

## Inputs for local functions:
type1 <- list(numeric = 1, numeric = 2, symm = 3, factor = 4)
type_name1 <- c("numeric", "numeric", "symm", "factor")
nvars <- length(type_name1)

## Format data nicely:
vars <- character()
cfs <- character()
sf <- ghs
data_type <- type1
data_name <- type_name1

for (z in 1:nvars) {
  var <- paste0("X", z)
  vars <- c(vars, var)
  cf <- paste0("CV", z)
  cfs[z] <- paste0(cf)
}

colnames(sf) <- c(vars, "X1_cat", "X2_cat")

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

ndom <- length(unique(sf$dom))

#------------------------ Register parallel backend --------------------------
n_cores <- 10

cl <- makeCluster(n_cores, outfile = "OUTPUT/parallel_log_bethel_atomic.txt")
registerDoParallel(cl)

#------------------------ Run Bethel allocation across seeds in parallel -------------
results <- foreach(
  seed = seeds,
  .packages = c("SamplingStrata", "dplyr"),
  .errorhandling = "pass"
) %dopar% {
  
  ## each worker needs its own copy of helper functions/env
  source("core_functions/medmix_gasampsi.R")
  source("core_functions/function_calc_variance.R")
  
  ########## LOAD INPUTS FROM MEDMIX (per seed)
  filename_in <- paste0(
    "OUTPUT/medmix_ga_atomic_", "seed", seed, "_results_LOCALTEST.Rdata"
  )
  load(filename_in)  # brings in `store_out` for this seed
  
  ssize <- sum(store_out$n)
  cv_y  <- unlist(store_out$cv_y_fpc)
  cv_b  <- unlist(store_out$cv_b_fpc)
  cv_g  <- unlist(store_out$cv_g_fpc)
  cv    <- cbind(t(cv_y), cv_b, cv_g)
  strata <- store_out$strata
  
  df_cv <- as.data.frame(cv)
  colnames(df_cv) <- c("cv_y_1", "cv_y_2", "cv_b_1", "cv_g_1")
  
  cv_list <- list(df_cv)
  error <- as.matrix(cv_list[[1]], nrow = 1)
  colnames(error) <- cfs
  
  ## specify error thresholds:
  error <- as.data.frame(list(
    DOM = rep("DOM1", ndom),
    error,
    domainvalue = c(1)
  ))
  
  ## ADD MEDMIX-SUGGESTED STRATA to this seed's own copy of the frame
  frame3_seed <- frame3
  frame3_seed$STRATO <- strata  # length must equal nrow(frame3), one label per unit
  
  strata_df <- frame3_seed %>%
    group_by(domainvalue, STRATO) %>%
    summarise(
      N  = n(),
      M1 = mean(Y1), S1 = sd(Y1),
      M2 = mean(Y2), S2 = sd(Y2),
      M3 = mean(Y3), S3 = sd(Y3),
      M4 = mean(Y4), S4 = sd(Y4),
      .groups = "drop"
    ) %>%
    mutate(
      COST = 1,
      CENS = 0,
      DOM1 = domainvalue
    ) %>%
    dplyr::select(-domainvalue)
  
  strata_df <- as.data.frame(strata_df)
  
  start.time <- Sys.time()
  proc.time.start <- proc.time()
  
  n <- tryCatch(
    bethel(
      stratif     = strata_df,
      errors      = error,
      minnumstrat = 2,
      printa      = TRUE
    ),
    error = function(e) {
      message(sprintf(
        "[seed=%d] bethel() failed: %s",
        seed, conditionMessage(e)
      ))
      stop(e)
    }
  )
  
  proc_time <- proc.time() - proc.time.start
  end.time <- Sys.time()
  time_taken <- difftime(end.time, start.time, units = "mins")
  
  n_strata_realized <- length(n)
  
  vars2 <- c("Y1", "Y2", "Y3", "Y4")
  fitness <- calculate_variance(
    df = frame3_seed, frame3_seed$STRATO,
    vars = vars2, n = n, type = data_type, type_name = data_name, seed = seed
  )
  
  # ----------------------- Store Results (per seed) --------------
  store <- list(
    max_clusters       = NA,
    seed               = seed,
    method             = "bethel (atomic)",
    nstrat_suggested   = store_out$nstrat,
    n_strata_realized  = n_strata_realized,
    df_sol             = frame3_seed,
    n                  = n,
    strata             = frame3_seed$STRATO,
    time               = time_taken,
    fitness            = fitness$deff,
    deff_y             = fitness$deff_y,
    deff_b             = fitness$deff_b,
    deff_g             = fitness$deff_g,
    cv_y               = fitness$cv_k,
    cv_b               = fitness$cv_b,
    cv_g               = fitness$cv_g,
    cv_y_fpc           = fitness$cv_k_fpc,
    cv_b_fpc           = fitness$cv_b_fpc,
    cv_g_fpc           = fitness$cv_g_fpc,
    proc_time          = proc_time,
    used_initial_solution = TRUE
  )
  
  filename <- paste0(
    "OUTPUT/bethel_ghs_atomic1_fpc_", "seed", seed, "_results_LOCALTEST.Rdata"
  )
  save(store, file = filename)
  
  store
}

stopCluster(cl)