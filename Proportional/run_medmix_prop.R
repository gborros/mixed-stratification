## Run PAM (parallel)

library(dplyr)
library(haven)
library(cluster)
library(doParallel)

#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

source("core_functions/medmix_gasampsi.R")
source("core_functions/function_calc_variance.R")
source("core_functions/function_sample_size.R")

#------------------------ Parallel setup ---------------------------------

ncores <- 40
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(cluster)
  source("core_functions/medmix_gasampsi.R")
  source("core_functions/function_calc_variance.R")
  source("core_functions/function_sample_size.R")
})

registerDoParallel(cl)

#------------------------ Section 1: Datasets ---------------------------------

load("multivar_datasets/d_mix1.RData")
d_mix1 <- as.data.frame(mix1)
type1 <- list(numeric=1, symm=2)
type_name1 <- c("numeric", "symm")

load("multivar_datasets/d_mix2.RData")
d_mix2 <- as.data.frame(mix2)
d_mix2$d_cat <- as.integer(d_mix2$d_cat == 2)
type2 <- list(numeric=1, symm=2, factor=3)
type_name2 <- c("numeric", "symm", "factor")

load("multivar_datasets/d_mix3.RData")
d_mix3 <- as.data.frame(mix3)
type3 <- list(numeric=1, symm=2)
type_name3 <- c("numeric", "symm")

load("multivar_datasets/d_mix4.RData")
d_mix4 <- as.data.frame(mix4)
d_mix4$d_cat <- as.integer(d_mix4$d_cat == 2)
type4 <- list(numeric=1, symm=2, factor=3)
type_name4 <- c("numeric", "symm", "factor")

load("multivar_datasets/d_mix5.RData")
d_mix5 <- as.data.frame(mix5)
type5 <- list(numeric=1, numeric=2, symm=3, symm=4)
type_name5 <- c("numeric", "numeric", "symm", "symm")

load("multivar_datasets/d_mix6.RData")
d_mix6 <- as.data.frame(mix6)
type6 <- list(numeric=1, numeric=2, symm=3, symm=4)
type_name6 <- c("numeric", "numeric", "symm", "symm")

datasets <- c("d_mix1", "d_mix2", "d_mix3", "d_mix4",
              "d_mix5", "d_mix6")

dfs <- list(d_mix1, d_mix2, d_mix3, d_mix4, d_mix5, d_mix6)

sample <- rep(750, 6)
for (num_strat in c(3, 4, 5, 6)) {

type <- list(type1, type2, type3, type4, type5, type6)
type_name <- list(type_name1, type_name2, type_name3,
                  type_name4, type_name5, type_name6)

#------------------------ Main loop ---------------------------------

for (dta in 1:length(datasets)) {
  
  ssize <- sample[dta]
  sf <- dfs[[dta]]
  data_type <- type[[dta]]
  data_name <- type_name[[dta]]
  num_strata <- num_strat
  
  vars <- paste0("X", 1:ncol(sf))
  colnames(sf) <- vars
  df <- sf
  
  #------------------------ Parallel iterations --------------------------
  
  results <- foreach(iter = 1:40,
                     .packages = c("dplyr", "cluster")) %dopar% {
                       
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
                       
                       list(
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

}

#------------------------ Stop cluster ---------------------------------

stopCluster(cl)