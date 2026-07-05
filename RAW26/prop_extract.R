## kmeans data extraction

library(dplyr)
setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed/RAW26/prop")

files <- list.files(
  "C:/Users/01459189/OneDrive/phd/hpc - mixed/RAW26/prop",
  pattern = "med_mixed_prop_.*\\.Rdata$",
  full.names = TRUE
)

parse_settings <- function(path) {
  name <- basename(path)
  
  dataset <- sub("^med_mixed_prop_(d_mix[0-9]+)_.*", "\\1", name)
  strata  <- as.numeric(sub(".*_([0-9]+)strata_.*", "\\1", name))
  N       <- as.numeric(sub(".*_([0-9]+)n_.*", "\\1", name))
  
  list(
    dataset = dataset,
    strata = strata,
    N = N
  )
}

res_list <- list()

base_out <- "C:/Users/01459189/OneDrive/phd/hpc - mixed/WIP26/prop/"

for (f in files) {
  
  load(f)  # loads object called `store`
  
  fname <- tools::file_path_sans_ext(basename(f))
  settings <- parse_settings(fname)
  
  # ---- extract cleanly ----
  n      <- store[[2]]
  strata <- store[[3]]
  deff     <- store[[4]]
  if (settings$dataset %in% c("d_mix5", "d_mix6")) {
    
    deff_y <- matrix(unlist(store[[5]]), ncol = 2, byrow = TRUE)
    deff_b <- matrix(unlist(store[[6]]), ncol = 2, byrow = TRUE)
    deff_g <- matrix(NA, nrow = 40, ncol = 2)
    
    cv_y <- matrix(unlist(store[[8]]), ncol = 2, byrow = TRUE)
    cv_b <- matrix(unlist(store[[9]]), ncol = 2, byrow = TRUE)
    cv_g <- matrix(NA, nrow = nsim, ncol = 2)
  }
  
  else {
    deff_y <- cbind(unlist(store[[5]]), NA)
    deff_b <- cbind(unlist(store[[6]]), NA)
    deff_g <- cbind(unlist(store[[7]]), NA)
    
    cv_y <- cbind(unlist(store[[8]]), NA)
    cv_b <- cbind(unlist(store[[9]]), NA)
    cv_g <- cbind(unlist(store[[10]]), NA)
  }
  
  if (ncol(deff_y) == 1) deff_y <- cbind(deff_y, NA)
  if (ncol(deff_b) == 1) deff_b <- cbind(deff_b, NA)
  if (ncol(deff_g) == 1) deff_g <- cbind(deff_g, NA)
  
  if (ncol(cv_y) == 1) deff_y <- cbind(cv_y, NA)
  if (ncol(cv_b) == 1) deff_b <- cbind(cv_b, NA)
  if (ncol(cv_g) == 1) deff_g <- cbind(cv_g, NA)
  
  time   <- store[[11]]

  # ---- save side objects ----
  save(n, file = paste0(base_out, fname, "_n.RData"))
  save(strata, file = paste0(base_out, fname, "_strata.RData"))
  save(deff_y, file = paste0(base_out, fname, "_deff_y.RData"))
  save(deff_b, file = paste0(base_out, fname, "_deff_b.RData"))
  save(deff_g, file = paste0(base_out, fname, "_deff_g.RData"))
  save(cv_y, file = paste0(base_out, fname, "_cv_y.RData"))
  save(cv_b, file = paste0(base_out, fname, "_cv_b.RData"))
  save(cv_g, file = paste0(base_out, fname, "_cv_g.RData"))
  
  # ---- clean tibble row ----
  res_tmp <- tibble(
    method  = "prop",
    dataset = settings$dataset,
    strata  = settings$strata,
    N       = settings$N,
    deff   = unname(as.vector(deff)),
    time    = unname(as.vector(time)),
    deff_y_1 = deff_y[,1],
    deff_y_2 = deff_y[,2],
    deff_b_1 = deff_b[,1],
    deff_b_2 = deff_b[,2],
    deff_g_1 = deff_g[,1],
    deff_g_2 = deff_g[,2],
    ## CV
    cv_y_1 = cv_y[,1],
    cv_y_2 = cv_y[,2],
    cv_b_1 = cv_b[,1],
    cv_b_2 = cv_b[,2],
    cv_g_1 = cv_g[,1],
    
    run_id  = fname
  )
  
  res_list[[fname]] <- res_tmp
  
  rm(store)  # avoid accidental reuse
}

# combine everything
res <- bind_rows(res_list)

save(res, file = paste0(base_out, "prop_all_res.RData"))