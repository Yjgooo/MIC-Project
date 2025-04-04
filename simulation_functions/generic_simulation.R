#separation of the sim_model and fit model
generic_sim <- function(n, m, setup, cov = NULL, sim_model = "simulate_data_np", model = "np") { 
  # Retrieve the simulation function based on the provided name
  sim <- get(sim_model)
  
  # Use a large sample to get the true distribution of T
  data_large <- sim(n = 10000, setup = setup, cov = cov)
  true_mean <- mean(data_large$T)
  
  full_support <- 1:10
  true_pmf <- prop.table(table(factor(data_large$T, levels = full_support)))
  
  vec <- numeric(m)
  tv_dist_vec <- numeric(m)
  
  # Models on restricted data
  for(i in 1:m) {
    # Note: cov parameter is not passed here to mimic the original function.
    data_oracle <- sim(n = n, setup = setup)
    data_observed <- subset(data_oracle, select = -T)  # dataframe without T
    df_ <- transform_df(data_observed, score_max = 11)
    
    # For new methods, do not pass user_formula
    if(model %in% c("roc", "pred", "pred.adj")){
      fit <- MIC(df_, method = model, cov = cov)
    } else {
      fit <- MIC(df_, user_formula = "cont + factor(sex) + factor(race)", method = model, cov = cov)
    }
    
    vec[i] <- abs(fit$meanT - true_mean)
    
    # For new methods, density estimates are not provided so TV distance is NA.
    if(model %in% c("roc", "pred", "pred.adj")){
      tv_dist_vec[i] <- NA
    } else {
      tv_dist_vec[i] <- tv_distance(true_pmf, pmf(fit$model, cov, max_time = 10)$y)
    }
  }
  
  result <- list(
    model = model, 
    dif = vec,
    tv_dist = tv_dist_vec,  # For new methods this will be NA
    n = n,
    m = m,
    setup = setup
  )
  
  class(result) <- "mega"
  return(result)
}

#setup is fixed 
run_simulations <- function(setup, cov = NULL, m = 200) {
  # Define the sample sizes and simulation model choices
  n_values <- c(100, 500, 1000)
  sim_models <- c("simulate_data_po", "simulate_data_ph", "simulate_data_np", "simulate_data_real")
  models <- c("np", "ph", "po", "pred", "pred.adj")
  
  # Container to store results for each combination
  results_list <- list()
  
  # Loop over each sample size and simulation/model choice
  for (n_val in n_values) {
    for (i in seq_along(sim_models)) {
      for (j in seq_along(models)) {
        sim_model <- sim_models[i]
        model <- models[j]
        
        # Call the generic simulation function
        sim_result <- generic_sim(n = n_val, m = m, setup = setup, cov = cov,
                                  sim_model = sim_model, model = model)
        
        # Create a key for the result list, e.g. "n_100_po"
        key <- paste("n", n_val, sim_model, "fit_model", model, sep = "_")
        results_list[[key]] <- sim_result 
      }
    }
  }
  
  return(results_list)
}

print.mega <- function(x) {
  cat(x$model, "Simulation Results\n")
  cat("------------------------\n")
  cat("Parameters: n =", x$n, ", m =", x$m, ", setup =", x$setup, "\n\n")
  
  cat("MSE:", mean((x$dif)^2, na.rm = TRUE), "\n")
  cat("95th percentile of MSE:", 
      as.numeric(quantile((x$dif)^2, probs = 0.95, na.rm = TRUE)), "\n")
  cat("MAD:", mean(abs(x$dif), na.rm = TRUE), "\n")
  cat("95th percentile of MADs:", 
      as.numeric(quantile(abs(x$dif), probs = 0.95, na.rm = TRUE)), "\n")
  
  if(all(is.na(x$tv_dist))) {
    cat("Average TV distance: Not available\n")
  } else {
    cat("Average TV distance:", mean(x$tv_dist, na.rm = TRUE), "\n")
  }
  
  invisible(x)
}