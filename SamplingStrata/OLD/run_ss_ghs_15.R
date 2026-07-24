## Run File

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)
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
num_strata = 20
ssize = 2000

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
error <- as.matrix(cv[[1]], nrow=1)
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
  
frame3 <- buildFrameDF(df = sf,
                         id = "id",
                         X = c("X1","X2", "X3", "X4"),
                         Y = c("X1","X2", "X3", "X4"),
                         domainvalue = "dom")
  
  
  
  ## specify error thresholds:
  ndom <- length(unique(sf$dom))
  error <- as.data.frame(list(DOM=rep("DOM1",ndom),
                              error,
                              domainvalue=c(1) ))

  
  #------------------------ Run GGA ----------------------------------
  store_n_all <- list()
  store_obj_all <- list()
  store_strata_all <- list()
  store_time_all <- list()
  store_fitness_all <- list()
  store_fitness_log_all <- list()
  
  results <- vector("list", 1)  # preallocate
  
  for (iter in 1:1) {
    set.seed(iter)
    start.time <- Sys.time()
    proc.time.start <- proc.time()
    solution <- optimStrata(method = "continuous",
                            errors = error, 
                            framesamp = frame3,
                            iter = 1000,
                            pops = 50,
                            nStrata = num_strata,
                            elitism_rate=0.2,
                            mut_chance = 0.7)
    proc_time  <- proc.time() - proc.time.start
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units="mins")
    
    strataStructure <- summaryStrata(solution$framenew,
                                     solution$aggr_strata,
                                     progress = FALSE)
    
    n <- as.vector(strataStructure$Allocation)
    df_sol <- solution$framenew
    df_sol$strata <- df_sol$STRATO
    
    
    fitness <- calculate_variance(df=df_sol, df_sol$strata, vars = vars, n=n, type=data_type, type_name=data_name, seed = iter)
    
    
    results[[iter]] <-     list(
      df_sol = df_sol,
      n = n,
      strata = df_sol$strata,
      time = time,
      fitness = fitness$deff,
      deff_y = fitness$deff_y,
      deff_b = fitness$deff_b,
      deff_g = fitness$deff_g,
      cv_y = fitness$cv_k,
      cv_b = fitness$cv_b,
      cv_g = fitness$cv_g,
      proc_time = proc_time
    )
  }
  
  store_df_sol_all      <- lapply(results, `[[`, "df_sol")
  store_n_all      <- do.call(rbind, lapply(results, `[[`, "n"))
  store_strata_all <- lapply(results, `[[`, "strata")
  store_time_all <- do.call(rbind, lapply(results, `[[`, "time"))
  store_fitness_all <- do.call(rbind, lapply(results, `[[`, "fitness"))
  store_deff_y_all <- lapply(results, `[[`, "deff_y")
  store_deff_b_all <- lapply(results, `[[`, "deff_b")
  store_deff_g_all <- lapply(results, `[[`, "deff_g")
  store_cv_y_all <- lapply(results, `[[`, "cv_y")
  store_cv_b_all <- lapply(results, `[[`, "cv_b")
  store_cv_g_all <- lapply(results, `[[`, "cv_g")
  store_proctime_all <- do.call(rbind, lapply(results, `[[`, "proc_time"))
  
  # ----------------------- Store Results --------------------------
  
  store <- list(
    store_df_sol_all,
    store_n_all,
    store_strata_all,
    store_time_all,
    store_fitness_all,
    store_deff_y_all,
    store_deff_b_all,
    store_deff_g_all,
    store_cv_y_all,
    store_cv_b_all,
    store_cv_g_all,
    store_proctime_all
  )
  
  filename = paste0("OUTPUT/ss_", datasets, "_", num_strata, "strata", "_", ssize, "n", "_results.Rdata")
  
  save(store, file = filename)