## Sample size

calculate_sample_size <- function(df,   ## data frame with y variables
                               strata, ## data frame with strata variable
                               vars,   ## vector of column names for y variables
                               method, ## allocation method
                               ssize,
                               type)  ## sample size
{
  df$strata <- strata
  dfg <- df %>%
    group_by(strata) %>%
    summarise(N_h = n(), .groups = "drop")
  
  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  k <- length(vars) ## number of variables

  if (method == 'prop') {
  
    n <- numeric()
    for (h in 1:L){
      n[h] <- as.integer(ssize*(N[h]/sum(N)))
    }
    
    for (h in 1:L){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    z <- which.max(N)
    n[z] <- n[z] + (ssize-sum(n)) ## allocate any remaining sample to largest stratum
  }
  
  if (method == 'gower') {
    
    V_strat <- calc_diss(df, type) ## calculate dissimilarity by stratum
    #V_tot <- sum(V_strat)
    #V_s_norm <- (V_strat/V_tot)*1
    
    n <- numeric()
    for (h in 1:L){
      #n[h] <- as.integer(ssize*(V_s_norm[h])) ## proportional to dissim
      n[h] <- as.integer(ssize*((N[h]*V_strat[h])/(sum(N*V_strat)))) ## neyman style
    }
    
    z <- which.max(V_strat)
    
    for (h in 1:L){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    n[z] <- n[z] + (ssize-sum(n)) ## allocating remaining sample to strata with most variance
    
  }
  
  
  return(n)
}

