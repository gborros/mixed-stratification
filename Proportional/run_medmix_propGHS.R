## Run PAM (parallel)

library(dplyr)
library(haven)
library(cluster)

#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

source("core_functions/medmix_gasampsi.R")
source("core_functions/function_calc_variance.R")
source("core_functions/function_sample_size.R")

#------------------------ Section 1: Datasets ---------------------------------

## GHS
load("multivar_datasets/ghs_2024_mixed.RData")
ghs <- as.data.frame(df)
ghs <- ghs[ghs$fin_reqinc!=9999999,]
ghs$fin_reqinc <- as.numeric(ghs$fin_reqinc)
ghs <- ghs[,-4]
ghs <- ghs[,-3]
ghs <- ghs[,-3]

ghs$hwl_status <- as.factor(ghs$hwl_status)
ghs$geotype <- as.factor(ghs$geotype)
ghs$head_age <- as.integer(ghs$head_age)

ghs$hwl_status <-as.integer(ghs$hwl_status == 5 | ghs$hwl_status == 6) ## wellbeing status 5 or 6
ghs$geotype <-as.integer(ghs$geotype==1) ## urban


type1 <- list(numeric=1, numeric=2, symm=3, factor=4)
type_name1 <- c("numeric", "numeric", "symm", "factor")

datasets <- c("ghs")
dfs <- list(ghs)
type <- list(type1)
type_name <- list(type_name1)
num_strata = 15
ssize = 2000

#------------------------ Main loop ---------------------------------

for (dta in 1:length(datasets)) {
  
  sf <- dfs[[dta]]
  data_type <- type[[dta]]
  data_name <- type_name[[dta]]
  
  vars <- paste0("X", 1:ncol(sf))
  colnames(sf) <- vars
  df <- sf
  
  results <- vector("list", 1)
  
  for (iter in 1:1) {
    
    set.seed(iter)
    start.time <- Sys.time()
    diss <- daisy(df, metric = "gower", type = data_type)
    result <- pam(diss, k = num_strata, diss = TRUE, nstart = 15)
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units = "mins")
    
    tmp <- data.frame(df, strata = result$clustering)
    strata <- as.factor(tmp$strata)
    
    n <- calculate_sample_size(
      df = tmp,
      strata = strata,
      vars = vars,
      method = "prop",
      ssize = ssize,
      type = data_type
    )
    
    calc <- calculate_variance(
      df = tmp,
      strata = strata,
      vars = vars,
      n = n,
      type = data_type,
      type_name = data_name,
      seed = iter
    )
    
    results[[iter]] <- list(   
                         n = n,
                         strata = strata,
                         deff = calc$deff,
                         deff_y = calc$deff_y,
                         deff_b = calc$deff_b,
                         deff_g = calc$deff_g,
			cv_k = calc$cv_k,
    			cv_b = calc$cv_b,
    			cv_g = calc$cv_g,
                         time = time
                       )
  }
  #------------------------ Combine results --------------------------
  
  store_n       <- do.call(rbind, lapply(results, `[[`, "n"))
  store_deff    <- do.call(rbind, lapply(results, `[[`, "deff"))
  store_deff_y  <- do.call(rbind, lapply(results, `[[`, "deff_y"))
  store_deff_b  <- do.call(rbind, lapply(results, `[[`, "deff_b"))
  store_deff_g  <- do.call(rbind, lapply(results, `[[`, "deff_g"))
  store_cv_y  <- do.call(rbind, lapply(results, `[[`, "cv_k"))
  store_cv_b  <- do.call(rbind, lapply(results, `[[`, "cv_b"))
  store_cv_g  <- do.call(rbind, lapply(results, `[[`, "cv_g"))
  store_time    <- do.call(rbind, lapply(results, `[[`, "time"))
  store_strata  <- lapply(results, `[[`, "strata")
  
  store <- list(df,
                store_n,
                store_strata,
                store_deff,
                store_deff_y,
                store_deff_b,
                store_deff_g,
                store_cv_y,
                store_cv_b,
                store_cv_g,
                store_time)
  
  filename <- paste0(
    "OUTPUT/",
    "med_mixed_prop_", datasets[dta], "_",
    num_strata, "strata_", ssize, "n_results.Rdata"
  )
  
  save(store, file = filename)
}
