## Calculate Variance Function

### Note: CV has been adjusted to include versions with and without FPC
### All other measures (variance, design effect) include the FPC.

calculate_variance <- function(df,   ## data frame with y variables
                               strata, ## data frame with strata variable
                               vars,   ## vector of column names for y variables
                               n,  ## sample allocation across strata 
                               type,
                               type_name, 
                               seed)      ## objective function (determinant of var-cov, total variance, sum of squares)
{
  df$strata <- strata
  dfg <- df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(N_h = n(), .groups = "drop")

  N <- dfg$N_h 
  W <- dfg$N_h / sum(N) ## population proportion of stratum
  L <- nrow(dfg) ## number of strata
  k <- length(vars) ## number of variables
  ssize <- sum(n)
  pop_tot <- sum(N)
  
  v_stk <- numeric()
  v_srsk <- numeric()
  deff_y <- numeric()
  v_stb <- numeric()
  v_srsb <- numeric()
  deff_b <- numeric()
  v_stg <- numeric()
  v_srsg <- numeric()
  deff_g <- numeric()
  cv_k <- numeric()
  cv_b <- numeric()
  cv_g <- numeric()
  cv_k_fpc <- numeric()
  cv_b_fpc <- numeric()
  cv_g_fpc <- numeric()
  
  for (j in 1:k) {
    
    ## MEAN
    if (type_name[j] == "numeric") {
      Y = df[[vars[j]]]
      
      ## SRS
      Y_m = as.matrix(mean(Y)) ## population mean
      V = stats::var(Y) ## population SD 
      
      #* Design based SRS variance:
      v_srs = V/ssize*((pop_tot-ssize)/pop_tot)
      
      df_st <- df %>%
        dplyr::mutate(Y_tmp = df[[vars[j]]]) %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
          N_h = dplyr::n(),
          var = ifelse(N_h == 1, 0, var(Y_tmp, na.rm = TRUE)),
          ybar = mean(Y_tmp, na.rm = TRUE),
          .groups = "drop"
        ) ## stratified SD (population); treating any singletons as 0 variance
      
      Vh <- df_st$var
      v_st <- numeric()
      
      #* Design based stratified sample variance:
      for (h in 1:L) {
        v_st[h] <- (W[h]^2 * (1 - (n[h] / N[h])) * (Vh[h] / n[h])) ## Kozak (5)
      }
      
      v_st <- sum(v_st)
      deff = v_st/v_srs
      
      #* Design based stratified sample mean:
      ybar <- df_st$ybar
      ybar_st <- sum(W * ybar)
      
      v_st_nofpc <- numeric()
      #* Design based stratified sample variance with NOFPC - for comparison to SamplingStrata:
      for (h in 1:L) {
        v_st_nofpc[h] <- (W[h]^2 * (Vh[h] / n[h])) ## Kozak (5)
      }
      
      v_st_nofpc <- sum(v_st_nofpc) ## sum over strata
      
      #* Design based stratified sample CV (NO FPC):
      cv_st = sqrt(v_st_nofpc)/ybar_st
      
      #* Design based stratified sample CV (WITH FPC):
      cv_st_fpc = sqrt(v_st)/ybar_st
      
      ## Collect output of interest:
      cv_k <- rbind(cv_k, cv_st)
      cv_k_fpc <- rbind(cv_k_fpc, cv_st_fpc)
      v_stk <- rbind(v_stk, v_st)
      v_srsk <- rbind(v_srsk, v_srs)
      deff_y = rbind(deff_y, deff)
    }
    
    ## PROPORTION
    if (type_name[j] == "symm") {
      B = df[[vars[j]]]
      
      ## SRS
      A = sum(B)
      P = A/pop_tot ## population proportion
      
      #* Design based SRS variance:
      v_srs = (P*(1-P)/ssize)*((pop_tot-ssize)/(pop_tot-1))
      
      ## Stratified sample
      df_st <- df %>%
        dplyr::mutate(B_tmp = df[[vars[j]]]) %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
        sum = sum(B_tmp, na.rm = TRUE),
          .groups = "drop"
        ) ## stratified SD
      
      
      Ph <- df_st$sum/N
      v_st <- numeric()
      
      #* Design based stratified sample variance:
      for (h in 1:L) {
        v_st[h] <- (1/pop_tot)^2*(N[h]^2*(N[h]-n[h])/(N[h]-1))*((Ph[h]*(1-Ph[h]))/n[h])
      }
      
      v_st <- sum(v_st)
      deff = v_st/v_srs
      
      #* Design based stratified sample proportion:
      phat_st = sum(W*Ph)
      
      v_st_nofpc <- numeric()
      #* Design based stratified sample variance with NO FPC (for comparison to SamplingStrata):
      for (h in 1:L) {
        v_st_nofpc[h] <- (1/pop_tot)^2 * N[h]^2 * ((Ph[h]*(1-Ph[h]))/n[h])
      }

      v_st_nofpc <- sum(v_st_nofpc)
      
      #* Design based stratified sample CV:
      cv_st <- sqrt(v_st_nofpc) / phat_st
      
      #* Design based stratified sample CV (WITH FPC):
      cv_st_fpc = sqrt(v_st)/phat_st
      
      # Collect output of interest:
      cv_b <- rbind(cv_b, cv_st)
      cv_b_fpc <- rbind(cv_b_fpc, cv_st_fpc)
      v_stb <- rbind(v_stb, v_st)
      v_srsb <- rbind(v_srsb, v_srs)
      deff_b = rbind(deff_b, deff)
    }
    
    ## TOTAL
    if (type_name[j] == "factor") {
      G = df[[vars[j]]]
      
      ## SRS
      A = sum(G)
      P = A/pop_tot ## population proportion
      
      #* Design based SRS variance:
      v_srs = pop_tot^2*(P*(1-P)/ssize)*((pop_tot-ssize)/(pop_tot-1))
      
      ## Stratified sample
      df_st <- df %>%
        dplyr::mutate(G_tmp = df[[vars[j]]]) %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
        sum = sum(G_tmp, na.rm = TRUE),
          .groups = "drop"
        ) ## stratified SD
      
      
      Ph <- df_st$sum/N
      v_st <- numeric()
      
      #* Design based stratified sample variance:
      for (h in 1:L) {
        v_st[h] <- (1/pop_tot)^2*(N[h]^2*(N[h]-n[h])/(N[h]-1))*((Ph[h]*(1-Ph[h]))/n[h])
      }
      
      v_st <- pop_tot^2*sum(v_st)
      deff = v_st/v_srs
      
      #* Design based stratified sample category total:
      that_st = sum(N*Ph)
      
      v_st_nofpc <- numeric()
      #* Design based stratified sample variance with NO FPC (for comparison to SamplingStrata):
      for (h in 1:L) {
        v_st_nofpc[h] <- (1/pop_tot)^2 * N[h]^2 * ((Ph[h]*(1-Ph[h]))/n[h])
      }
      
      v_st_nofpc <- pop_tot^2*sum(v_st_nofpc)
      
      #* Design based stratified sample category CV:
      cv_st = sqrt(v_st_nofpc)/that_st
      
      #* Design based stratified sample CV (WITH FPC):
      cv_st_fpc = sqrt(v_st)/that_st
      
      # Collect output of interest:
      cv_g = rbind(cv_g, cv_st)
      cv_g_fpc <- rbind(cv_g_fpc, cv_st_fpc)
      v_stg <- rbind(v_stg, v_st)
      v_srsg <- rbind(v_srsg, v_srs)
      deff_g = rbind(deff_g, deff)
    }
      
    }
 
  deff_tot <- sum(deff_y) + sum(deff_b) + sum(deff_g)
  
  return(list(
    cv_k = cv_k,
    cv_b = cv_b,
    cv_g = cv_g,
    cv_k_fpc = cv_k_fpc,
    cv_b_fpc = cv_b_fpc,
    cv_g_fpc = cv_g_fpc,
    deff_y   = deff_y,
    deff_b   = deff_b,
    deff_g   = deff_g,
    deff = deff_tot
  ))
}   
  
  
