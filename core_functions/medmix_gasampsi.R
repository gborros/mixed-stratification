## Clustering with GA sample size 

############################################################
############################################################

# ---- Section 2: Supporting Functions for GA ----

# Initialise population using mclust and random sample size allocation
init_pop_gga <- function(pop_size, num_strata, n_obs, ssize, seed, df, type, type_name) {
  population <- list()
  sampsi <- list()
  
  set.seed(seed)
  ## Generate boundaries once, afterwhich the sample will be randomly generated
  diss <- daisy(df, metric = "gower", type=type)
  result <- pam(diss, k = num_strata, diss = TRUE, nstart=15)
  df <- data.frame(df, result$clustering)
  try <- 1
  while (length(population) < pop_size) {
    ##########* Stratification *##########
    ## Initial population based on clustering
    try = try + 1
    set.seed(seed+length(population)+try)
    candidate <- as.numeric(result$clustering)
    
    ##########* Sample Allocation *########
    length <- numeric()
    for (h in 1:num_strata){
      length <- rbind(length, length(candidate[candidate==h]))
    }
    
    ## Random sample allocation for initialisation into GGA
    n <- numeric()
    samp_prop <- floor(runif(num_strata,1,10)) 
    sum_samp <- sum(samp_prop)
    n <- floor(samp_prop/sum_samp*ssize)
    
    if (sum(n)>ssize) {
      N=sum(n)
      n[which.min(N-n)] <- n[which.min(N-n)] - (sum(n)-ssize) 
    }
    
    if (sum(n)<ssize) {
      N=sum(n)
      n[which.max(N-n)] <- n[which.max(N-n)] + (ssize-sum(n)) 
    }
    
    for (h in 1:num_strata){
      if (n[h] %in% c(0,1)) { ## ensure smallest stratum is bigger than 1 
        n[h] = 2
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    for (h in 1:num_strata){
      if (n[h] > length[h]) { ## ensure sample is smaller than pop
        n[h] = length[h]
        n[which.max(n)] <- ssize-sum(n)+n[which.max(n)]
      }
    }
    
    N=sum(n)
    z <- which.max(N)
    n[z] <- n[z] + (ssize-sum(n)) ## allocate any remaining sample to largest stratum
    
    ## Constraint check
    c1 <- c()
    for (h in 1:num_strata){
      c1[h] <- n[h] <= length[h] & 2 <= n[h]}
    
    n <- n
    
    # Check groups: must have two observations per group and cover all groups
    group_counts <- table(candidate)
    all_groups_present <- length(group_counts) == num_strata
    all_groups_have_min_two <- all(group_counts >= 2)
    all_c1_true <- all(c1)==TRUE
    all_groups_meet_ssize <- sum(n) == ssize
    
    # If all groups have at least 2 observations, accept the candidate
    if (all_groups_present && all_groups_have_min_two && all_c1_true && all_groups_meet_ssize) {
      sampsi[[length(population) + 1]] <- n
      population[[length(population) + 1]] <- candidate
    } else {
      # If any group has less than 2 observations, redo the sample
      print("Retrying initial sample")
      next
    }
  }
  
  
  return(list(population=population, 
              n_list=sampsi))
}

#init <- init_pop_gga(25, 3, 146, 50)
#init <- init_pop_gga(pop_size=25, num_strata=3, n_obs=5000, ssize=500, seed=1, df=df, type=data_type, type_name=data_name)
############################################################

evaluate_fitness_gga <- function(df, population, vars, ssize, type, type_name, seed) {
  set.seed(seed)
  n_list <- population$n_list
  population <- population$population
  
  # Step 1: keep full results
  results <- mapply(function(individual, n) {
    calculate_variance(df, individual, vars, n, type, type_name, seed)
  }, population, n_list, SIMPLIFY = FALSE)
  
  # Step 2: extract components
  fitness <- sapply(results, function(x) x$deff)
  deff_y  <- lapply(results, function(x) x$deff_y)
  deff_b  <- lapply(results, function(x) x$deff_b)
  deff_g  <- lapply(results, function(x) x$deff_g)
  cv_y  <- lapply(results, function(x) x$cv_k)
  cv_b  <- lapply(results, function(x) x$cv_b)
  cv_g  <- lapply(results, function(x) x$cv_g)
  
  return(list(
    fitness = -fitness,
    n       = n_list,
    deff_y  = deff_y,
    deff_b  = deff_b,
    deff_g  = deff_g,
    cv_y = cv_y,
    cv_b = cv_b,
    cv_g = cv_g
  ))
}


#fit <- evaluate_fitness_gga(df=df, population=init, vars=vars, ssize=500, type=data_type, type_name=data_name, seed=1)


# Selection based on fitness
select_parents_gga <- function(population, fitness, num_parents) {
  
  n_parent <- population$n_list
  population <- population$population
  
  # Rank fitness values - trying hard ranking now, commenting out probability ranking
  #ranks <- rank(fitness, ties.method = "first") ## first ranked is the worst, last rank is best
  #prob <- ranks / sum(ranks)  # Calculate selection probabilities based on ranks
  #print(prob)
  # Sample parents based on the selection probabilities
  #selected_indices <- sample(1:length(population), size = num_parents, prob = prob, replace = FALSE)
  #selected_indices <- ranks[1:num_parents]
  selected_indices <- order(fitness, decreasing = TRUE)[1:num_parents]
  #print(selected_indices)
  # Return the selected parents from the valid population
  return(list(
    
    population[selected_indices],
    
    fitness[selected_indices],
    
    n_parent[selected_indices]
  ))
}

#calc_parents <- select_parents_gga(init, fit$fitness, 3)
#parents <- calc_parents[[1]]
#parent_fitness <- calc_parents[[2]]
#parents_n <- calc_parents[[3]]

crossover_gga <- function(parents, num_strata, parents_n, ssize, crossover_rate, seed, generation) {
  offspring <- list()
  offspring_n <- list()
  
  for (i in seq(1, length(parents) - 1, by = 2)) {
    p1 <- parents[[i]]
    p2 <- parents[[i + 1]]
    n1 <- parents_n[[i]]
    n2 <- parents_n[[i+1]]
    length <- numeric()
    set.seed(seed+i+generation)
    crossover_prob <- runif(1)
    #print(crossover_prob)
    
    if (crossover_prob >= crossover_rate) {
      
      child1 <- p1
      child_n1 <- n1
      child2 <- p2
      child_n2 <- n2
      
      offspring <- c(offspring, list(child1, child2))
      offspring_n <- c(offspring_n, list(child_n1, child_n2))
    }
    
    else {
      # -------------------------------- Creating Child 1 ------------------------
      groups <- unique(p1)
      selected_groups <- sample(groups, size = ceiling(length(groups) / 2))
      
      child1 <- rep(NA, length(p1))
      child_n1 <- rep(NA, length(n1))
      
      ################### Group assignment mutation #############################
      
      #child1[p1 %in% selected_groups] <- p1[p1 %in% selected_groups]
      #child1[is.na(child1)] <- p2[is.na(child1)]
      
      ################### Sample size per stratum mutation #######################
      
      child_n1[selected_groups] <- n1[selected_groups]
      child_n1[is.na(child_n1)] <- n2[is.na(child_n1)]
      
      ################### Total sample size checks and adjustments ###############
      
      if (sum(child_n1)>ssize) {
        #print("Sample size for child 1 is greater than required, performing adjustment")
        N=sum(child_n1)
        child_n1[which.min(N-child_n1)] <- child_n1[which.min(N-child_n1)] - (sum(child_n1)-ssize) 
      }
      
      if (sum(child_n1)<ssize) {
        #print("Sample size for child 1 is smaller than required, performing adjustment")
        N=sum(child_n1)
        child_n1[which.max(N-child_n1)] <- child_n1[which.max(N-child_n1)] + (ssize-sum(child_n1)) 
      }
      
      #print("Sample size for first child")
      #print(sum(child_n1))
      
      ################### Sample size per stratum checks #########################
      # Length (N_h) check for sample size
      for (h in 1:num_strata){
        length <- rbind(length, length(p1[p1==h]))
      }
      
      ## Constraint check
      c1 <- c()
      for (h in 1:num_strata){
        c1[h] <- child_n1[h] <= length[h] & 2 <= child_n1[h]}
      
      ## sample sizes that don't work should be unchanged:
      if(all(c1)==TRUE) {child_n1<-child_n1} else {child_n1<-n1} 
      #if(all(c1)==TRUE) {print("Sample size per stratum constraints met for child 1")} else {print("Sample size per stratum constraints not met for child 1, no crossover applied.")} 
      ###########################################################################
      
      group_counts <- table(p1)
      all_groups_present <- length(group_counts) == num_strata
      all_groups_have_min_two <- all(group_counts >= 2)
      all_groups_meet_ssize <- sum(child_n1) == ssize
      
      # If all groups have at least 2 observations, accept the child
      if (isTRUE(all_groups_present) && isTRUE(all_groups_have_min_two) && isTRUE(all_groups_meet_ssize)) {
        child1 <- p1
        child_n1 <- child_n1
      } else {
        # If child1 does not meet conditions, then leave unchanged
        child1 <- p1
        child_n1 <- n1
      }
      
      # ------------------------------- Creating Child 2 ------------------------
      groups2 <- unique(p2)
      selected_groups2 <- sample(groups2, size = ceiling(length(groups2) / 2))
      
      ################### Group assignment mutation #############################
      
      child2 <- rep(NA, length(p2))
      # child2[p2 %in% selected_groups2] <- p2[p2 %in% selected_groups2]
      # child2[is.na(child2)] <- p1[is.na(child2)]
      
      ################### Sample size per stratum mutation #######################
      
      child_n2 <- rep(NA, length(n2))
      
      child_n2[selected_groups2] <- n2[selected_groups2]
      child_n2[is.na(child_n2)] <- n1[is.na(child_n2)]
      
      ################### Total sample size checks and adjustments ###############
      
      if (sum(child_n2)>ssize) {
        # print("Sample size for child 2 is greater than required, performing adjustment")
        N=sum(child_n2)
        child_n2[which.min(N-child_n2)] <- child_n2[which.min(N-child_n2)] - (sum(child_n2)-ssize) 
      }
      
      if (sum(child_n2)<ssize) {
        #  print("Sample size for child 2 is smaller than required, performing adjustment")
        N=sum(child_n2)
        child_n2[which.max(N-child_n2)] <- child_n2[which.max(N-child_n2)] + (ssize-sum(child_n2)) 
      }
      
      # print("Sample size for second child")
      # print(sum(child_n2))
      
      ################### Sample size per stratum checks #########################
      # Length (N_h) check for sample size
      length <- numeric()
      for (h in 1:num_strata){
        length <- rbind(length, length(p2[p2==h]))
      }
      
      ## Constraint check
      c1 <- c()
      for (h in 1:num_strata){
        c1[h] <- child_n2[h] <= length[h] & 2 <= child_n2[h]}
      
      ## sample sizes that don't work should be unchanged:
      if(all(c1)==TRUE) {child_n2<-child_n2} else {child_n2<-n2} 
      #if(all(c1)==TRUE) {print("Sample size per stratum constraints met for child 2")} else {print("Sample size per stratum constraints not met for child 2, no crossover applied.")} 
      ###########################################################################
      
      group_counts <- table(p2)
      all_groups_present <- length(group_counts) == num_strata
      all_groups_have_min_two <- all(group_counts >= 2)
      all_groups_meet_ssize <- sum(child_n2) == ssize
      
      # If all groups have at least 2 observations, accept the child
      if (isTRUE(all_groups_present) && isTRUE(all_groups_have_min_two) && isTRUE(all_groups_meet_ssize)) {
        child2 <- p2
        child_n2 <- child_n2
      } else {
        # If child2 does not meet conditions, then leave unchanged
        child2 <- p2
        child_n2 <- n2
      }
      offspring <- c(offspring, list(child1, child2))
      offspring_n <- c(offspring_n, list(child_n1, child_n2))
    }
  }
  # Handle odd number of parents
  if (length(parents) %% 2 == 1) {
    offspring <- c(offspring, list(parents[[length(parents)]]))
    offspring_n <- c(offspring_n, list(parents_n[[length(parents_n)]]))
  }
  return(list(offspring=offspring,
              offspring_n=offspring_n))
}

#offspring <- crossover_gga(parents, 3, parents_n, 50, crossover_rate=0.9)

# Mutation by randomly changing group assignment 
mutate_population_gga <- function(population, mutation_rate, num_strata, offspring_n, seed, generation) {
  mutated_population <- list()
  mutated_n_list <- list()
  
  for (i in seq_along(population)) {
    individual <- population[[i]]
    n <- offspring_n[[i]]
    length <- numeric()
    seedling <- i+generation
    set.seed(seed+(seedling))
    if (runif(1) < mutation_rate) {
      
      # # Group assignment mutation
      # for (j in seq_along(individual)) {
      #   if (runif(1) < mutation_rate) {
      #     individual[j] <- sample(1:num_strata, 1)
      #   }
      # }
      
      # Length (N_h) calc for sample size
      for (h in 1:num_strata){
        length <- rbind(length, length(individual[individual==h]))
      }
      
      # ssize mutation
      adjust <- max(1, min(n)/2)  # don't want to have n[h]<2 or larger than N[h]
      X <- sample(1:adjust, 1)
      strata_indices <- sample(1:num_strata, 2, replace = FALSE)
      
      increase_idx <- strata_indices[1]
      decrease_idx <- strata_indices[2]
      
      n[decrease_idx] <- n[decrease_idx] - X
      n[increase_idx] <- n[increase_idx] + X
      
      ## Constraint check
      c1 <- c()
      for (h in 1:num_strata){
        c1[h] <- n[h] <= length[h] & 2 <= n[h]}
      
      ## sample sizes that don't work should be unchanged:
      if(all(c1)==TRUE) {n<-n} else {n<-offspring_n[[i]]} 
      
    }
    # Save mutated individual and sample size
    mutated_population[[i]] <- individual
    mutated_n_list[[i]] <- n
  }
  return(list(
    population = mutated_population,
    n_list = mutated_n_list
  ))
}


#mut <- mutate_population_gga(offspring$offspring, 0.5, 3, offspring$offspring_n); mut

fitness_log <- data.frame(
  Generation = integer(),
  BestFitness = numeric(),
  MeanFitness = numeric()
)


run_medmix_gasampsi <- function(seed, df, pop_size, num_strata, num_generations, mutation_rate, vars, ssize, elitism_rate, crossover_rate, type, type_name) {
  ## think about elitism
  n_obs <- nrow(df)
  population <- init_pop_gga(pop_size, num_strata, n_obs, ssize, seed, df, type, type_name)
  
  for (generation in 1:num_generations) {
    # Evaluate fitness
    fit <- evaluate_fitness_gga(df, population, vars, ssize, type, type_name, seed)
    fitness <- fit$fitness
    gen_n <- fit$n
    
    # Sort population by fitness (best first) and select elite according to elitism_rate
    sorted_idx <- order(fitness, decreasing = TRUE)
    num_elite <- floor(elitism_rate * pop_size)
    pop <- population$population
    elite_population <- pop[sorted_idx[1:num_elite]]
    elite_n <- gen_n[sorted_idx[1:num_elite]]
    
    # Selection: Select based on fitness for crossover later
    calc_parents <- select_parents_gga(population, fitness, pop_size - num_elite)
    parents <- calc_parents[[1]]
    parents_fitness <- calc_parents[[2]]
    parents_n <- calc_parents[[3]]
    
    # Pass parents and their fitness to crossover
    offspring <- crossover_gga(parents, num_strata, parents_n, ssize, crossover_rate, seed, generation)
    offspring_n <- offspring$offspring_n
    offspring<-offspring$offspring
    
    # Mutation
    mutation <- mutate_population_gga(offspring, mutation_rate, num_strata, offspring_n, seed, generation)
    offspring <- mutation$population
    offspring_n <- mutation$n_list
    
    # Combine elite solutions with new offspring
    population <- c(elite_population, offspring)
    new_n <- c(elite_n, offspring_n)
    
    population <- list(population=population,
                       n_list=new_n)
    
    ## Evaluate fitness for new pop
    fit <- evaluate_fitness_gga(df, population, vars, ssize, type, type_name, seed)
    fitness <- fit$fitness
    n <- fit$n
    best_idx <- which.max(fitness)
    best_fitness <- fitness[best_idx]
    best_n <- n[[best_idx]]
    
    #Print current best fitness
    cat("Generation", generation, "Best Fitness:", best_fitness, "\n")
    #Print current best sample size
    cat("Generation", generation, "Best Sample Allocation:", best_n, "\n")
    
    # Compute stats for graph
    best_fitness <- -max(fitness)
    mean_fitness <- -mean(fitness)
    
    # Log stats
    fitness_log <- rbind(fitness_log, data.frame(
      Generation = generation,
      BestFitness = best_fitness,
      MeanFitness = mean_fitness
    ))
    
    # # Live ggplot update
    # p <- ggplot(fitness_log, aes(x = Generation)) +
    #   geom_line(aes(y = BestFitness, color = "Best Fitness"), size = 1) +
    #   geom_line(aes(y = MeanFitness, color = "Mean Fitness"), linetype = "solid", size = 0.5) +
    #   labs(title = "GGA Fitness Progression", y = "Fitness") +
    #   scale_color_manual(values = c("Best Fitness" = "blue", "Mean Fitness" = "firebrick2")) +
    #   theme_minimal()
    # 
    # print(p)
    # Sys.sleep(0.1)  # brief pause so the plot window updates
    
  }
  
  # Final evaluation
  fit <- evaluate_fitness_gga(df, population, vars, ssize, type, type_name, seed)
  fitness <- fit$fitness
  best_idx <- which.max(fitness)
  n <- fit$n
  population <- population$population
  deff_y <- fit$deff_y
  deff_b <- fit$deff_b
  deff_g <- fit$deff_g
  cv_y <- fit$cv_y
  cv_b <- fit$cv_b
  cv_g <- fit$cv_g
  return(list(best_solution = population[[best_idx]], best_fitness = fitness[best_idx], best_n = n[[best_idx]], fitness_log = fitness_log, best_deff_y = deff_y[best_idx], best_deff_b = deff_b[best_idx], best_deff_g = deff_g[best_idx],
              best_cv_y = cv_y[best_idx], best_cv_b = cv_b[best_idx], best_cv_g = cv_g[best_idx]))
}

