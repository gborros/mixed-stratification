## Run File - SINGLE (num_strata, seed) COMBINATION, SEQUENTIAL, VERBOSE
## Continuous-only variant (2 numeric variables: X1, X2)
## Local testing version - no parallelism, max visibility into each step.
## Derived from the parallel HPC run file; same logic, same function calls,
## just run once with explicit prints so you can inspect intermediate objects.

library(dplyr)
library(haven)
library(cluster)
library(SamplingStrata)

#setwd("C:/Users/01459189/OneDrive/phd/hpc - mixed")

# ------------------------ SET YOUR TEST COMBINATION HERE -------------------
num_strata <- 6     # <-- change this
seed       <- 5      # <-- change this
ssize      <- 2000

cat("========================================\n")
cat(sprintf("RUN: num_strata = %d, seed = %d, ssize = %d\n", num_strata, seed, ssize))
cat("========================================\n\n")

# ------------------------ Source local functions ---------------------------
source("core_functions/medmix_gasampsi.R")
source("core_functions/function_calc_variance.R")

#------------------------ Data ---------------------------------
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
ghs <- ghs[, -4]
ghs <- ghs[, -3]
cat(sprintf("   After column drops: %d cols -> %s\n", ncol(ghs), paste(colnames(ghs), collapse = ", ")))

ghs$head_age <- as.integer(ghs$head_age)

cat("   Post-recode summary:\n")
print(summary(ghs))

type1 <- list(numeric = 1, numeric = 2)
type_name1 <- c("numeric", "numeric")

## CVs
cat("\n>> Loading CV input table...\n")
load("CVs/cv_input_app.RData")

df_cv <- df_cv_min %>%
  filter(
    type == "original",
    strata == num_strata,
    N == ssize
  ) %>%
  select(cv_y_1, cv_y_2)

cat(sprintf("   Matched %d row(s) of df_cv for strata=%d, N=%d:\n", nrow(df_cv), num_strata, ssize))
print(df_cv)

if (nrow(df_cv) == 0) {
  stop(sprintf("No CV row found for num_strata=%d, ssize=%d — check cv_input_app.RData filter.", num_strata, ssize))
}

vars <- character()
cfs <- character()
sf <- ghs
cv <- list(df_cv)
error <- as.matrix(cv[[1]], nrow = 1)
data_type <- type1
data_name <- type_name1

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

cat("\n>> Building sampling frame with buildFrameDF()...\n")
frame3 <- buildFrameDF(
  df = sf,
  id = "id",
  X = c("X1", "X2"),
  Y = c("X1", "X2"),
  domainvalue = "dom"
)
cat(sprintf("   Frame built: %d rows x %d cols\n", nrow(frame3), ncol(frame3)))
str(frame3)

## specify error thresholds:
ndom <- length(unique(sf$dom))
error <- as.data.frame(list(
  DOM = rep("DOM1", ndom),
  error,
  domainvalue = c(1)
))
cat("\n>> Error/CV threshold table passed to optimStrata:\n")
print(error)

#------------------------ Kmeans initial solution (optional) ------------
cat("\n>> Attempting KmeansSolution2() initial solution...\n")
set.seed(seed)
initial_solution <- tryCatch(
  {
    sol <- KmeansSolution2(
      frame = frame3,
      errors = error,
      nstrat = num_strata
    )

    n_groups <- length(unique(sol))
    cat(sprintf("   KmeansSolution2 returned %d distinct strata (expected %d).\n", n_groups, num_strata))

    if (n_groups != num_strata) {
      message(sprintf(
        "[num_strata=%d, seed=%d] KmeansSolution2 produced %d distinct strata (expected %d) — discarding suggestion.",
        num_strata, seed, n_groups, num_strata
      ))
      NULL
    } else {
      cat("   Kmeans initial solution:\n")
      print(sol)
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

#------------------------ Run GGA (single seed, single run) -----------------
cat("\n>> Running optimStrata()...\n")
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
  mut_chance = 0.3
)

if (!is.null(initial_solution)) {
  optimStrata_args$suggestions <- initial_solution
}

used_initial_flag <- !is.null(initial_solution)
cat(sprintf("   used_initial_solution = %s\n", used_initial_flag))

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
cat(sprintf("\n>> optimStrata finished in %.2f minutes.\n", as.numeric(time_taken)))

cat("\n>> Structure of solution object:\n")
str(solution, max.level = 1)

strataStructure <- summaryStrata(solution$framenew,
                                  solution$aggr_strata,
                                  progress = FALSE
)
cat("\n>> summaryStrata() output:\n")
print(strataStructure)

n <- as.vector(strataStructure$Allocation)
n_strata_realized <- length(n) # actual number of strata this run produced
cat(sprintf("\n>> Requested nStrata = %d | Realized strata = %d\n", num_strata, n_strata_realized))
cat(">> Allocation vector (n):\n")
print(n)
cat(sprintf(">> Total allocated sample size: %d\n", sum(n)))

df_sol <- solution$framenew
df_sol$strata <- df_sol$STRATO

cat("\n>> Strata frequency table (from df_sol$strata):\n")
print(table(df_sol$strata))

cat("\n>> Calculating variance / CV diagnostics via calculate_variance()...\n")
fitness <- calculate_variance(
  df = df_sol, df_sol$strata,
  vars = vars, n = n, type = data_type, type_name = data_name, seed = seed
)

cat("\n>> Fitness / CV results:\n")
cat(sprintf("   deff (overall) : %s\n", fitness$deff))
cat(sprintf("   deff_y         : %s\n", paste(fitness$deff_y, collapse = ", ")))
cat(sprintf("   deff_b         : %s\n", paste(fitness$deff_b, collapse = ", ")))
cat(sprintf("   deff_g         : %s\n", paste(fitness$deff_g, collapse = ", ")))
cat(sprintf("   cv_y (cv_k)    : %s\n", paste(fitness$cv_k, collapse = ", ")))
cat(sprintf("   cv_b           : %s\n", paste(fitness$cv_b, collapse = ", ")))
cat(sprintf("   cv_g           : %s\n", paste(fitness$cv_g, collapse = ", ")))

# ----------------------- Store Results (single run) --------------
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
  "OUTPUT/ss_ghs_cont_", num_strata, "strata_seed", seed, "_", ssize, "n", "_results_LOCALTEST.Rdata"
)

# Uncomment to actually write output during local testing:
# save(store, file = filename)
cat(sprintf("\n>> (Not auto-saving) Would have written: %s\n", filename))
cat(">> Uncomment the save() line above if you want the .Rdata written to disk.\n")

cat("\n========================================\n")
cat("DONE.\n")
cat("========================================\n")
