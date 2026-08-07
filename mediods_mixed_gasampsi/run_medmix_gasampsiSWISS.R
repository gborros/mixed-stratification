## Run File

library(dplyr)
library(haven)
library(cluster)
library(doParallel)
library(SamplingStrata)
setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")
source("core_functions/medmix_gasampsi.R") ## sourcing GGA function
source("core_functions/function_calc_variance.R")

ncores <- 40
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(dplyr)
  library(haven)
  library(cluster)
  source("core_functions/medmix_gasampsi.R")
  source("core_functions/function_calc_variance.R")
})
registerDoParallel(cl)

#------------------------ Data ---------------------------------

## SWISS

data(swissmunicipalities)
swissmun <- swissmunicipalities[swissmunicipalities$REG < 4,
                                c("Airbat","POPTOT")]

frame1 <- buildFrameDF(df = swissmun,
                       id = "COM",
                       X = c("POPTOT.cat","HApoly.cat"),
                       Y = c("Airbat","Surfacesbois"),
                       domainvalue = "REG")


num_strata = 27
ssize = 299
iter = 1

df = swissmun

result <- run_medmix_gasampsi(seed=iter, df = df, pop_size = 50, num_strata = num_strata,
                                  num_generations = 2000, mutation_rate = 0.7, 
                                  vars=vars, ssize=ssize, elitism_rate = 0.2,
                                  crossover_rate=0.7, type=data_type, type_name=data_name)
    proc_time  <- proc.time() - proc.time.start
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units="mins")
    # ----------------------- Report Results --------------------------
    
    strata <- result$best_solution
    df$strata <- strata
    
    list(
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
}
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

stopCluster(cl)