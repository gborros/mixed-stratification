## Run File - GHS with no skewed variable, normal only

library(dplyr)
library(haven)
library(cluster)
#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")
source("core_functions/medmix_gasampsi.R") ## sourcing GGA function
source("core_functions/function_calc_variance.R")

#------------------------ Data ---------------------------------

## GHS
load("multivar_datasets/ghs_2024_mixed.RData")
ghs <- as.data.frame(df)
ghs <- ghs[ghs$fin_reqinc!=9999999,]
ghs$fin_reqinc <- as.numeric(ghs$fin_reqinc)
ghs <- ghs[,-4]
ghs <- ghs[,-3]
ghs <- ghs[,-3]
ghs <- ghs[,-2]

ghs$hwl_status <- as.factor(ghs$hwl_status)
ghs$geotype <- as.factor(ghs$geotype)
ghs$head_age <- as.integer(ghs$head_age)

ghs$hwl_status <-as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
ghs$geotype <-as.integer(ghs$geotype==1) ## urban


type1 <- list(numeric=1, symm=2, factor=3)
type_name1 <- c("numeric", "symm", "factor")

datasets <- c("ghs_norm")
dfs <- list(ghs)
type <- list(type1)
type_name <- list(type_name1)
num_strata = 20
ssize = 1000

for (dta in 1:length(datasets)) {
  vars <- character()
  sf <- dfs[[dta]]
  data_type <- type[[dta]]
  data_name <- type_name[[dta]]
  
  for (z in 1:ncol(sf)) {
    
    var <- paste0("X", z)
    vars <- c(vars, var)
    
  }
  
  colnames(sf) <- vars
  df <- sf
  
  #------------------------ Run GGA ----------------------------------
  store_n_all <- list()
  store_obj_all <- list()
  store_strata_all <- list()
  store_time_all <- list()
  store_fitness_all <- list()
  store_fitness_log_all <- list()
  
  results <- vector("list", 1)  # preallocate
  
  for (iter in 1:1) {
    start.time <- Sys.time()
    proc.time.start <- proc.time()
    result <- run_medmix_gasampsi(
      seed = iter, df = df, pop_size = 50, num_strata = num_strata,
      num_generations = 1000, mutation_rate = 0.7, 
      vars = vars, ssize = ssize, elitism_rate = 0.2,
      crossover_rate = 0.7, type = data_type, type_name = data_name
    )
    
    proc_time  <- proc.time() - proc.time.start
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units="mins")
    
    strata <- result$best_solution
    
    results[[iter]] <-     list(
      n = result$best_n,
      strata = strata,
      time = time,
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
  }
  
  store_n_all      <- do.call(rbind, lapply(results, `[[`, "n"))
  store_strata_all <- lapply(results, `[[`, "strata")
  store_time_all <- do.call(rbind, lapply(results, `[[`, "time"))
  store_fitness_all <- do.call(rbind, lapply(results, `[[`, "fitness"))
  store_fitness_log_all <- lapply(results, `[[`, "fitness_log")
  store_deff_y_all <- lapply(results, `[[`, "deff_y")
  store_deff_b_all <- lapply(results, `[[`, "deff_b")
  store_deff_g_all <- lapply(results, `[[`, "deff_g")
  store_cv_y_all <- lapply(results, `[[`, "cv_y")
  store_cv_b_all <- lapply(results, `[[`, "cv_b")
  store_cv_g_all <- lapply(results, `[[`, "cv_g")
  store_proctime_all <- do.call(rbind, lapply(results, `[[`, "proc_time"))

# ----------------------- Store Results --------------------------

store <- list(
  df,
  store_n_all,
  store_strata_all,
  store_time_all,
  store_fitness_all,
  store_fitness_log_all,
  store_deff_y_all,
  store_deff_b_all,
  store_deff_g_all,
  store_cv_y_all,
  store_cv_b_all,
  store_cv_g_all,
  store_proctime_all
)

filename = paste0("OUTPUT/medmix_ga_", datasets[dta], "_", num_strata, "strata", "_", ssize, "n", "_results.Rdata")

save(store, file = filename)

}