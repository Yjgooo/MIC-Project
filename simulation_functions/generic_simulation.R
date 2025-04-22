generic_sim <- function(n, m, setup, cov = NULL,
                        sim_model = "simulate_data_np",
                        model     = "np",
                        N0        = 10000,
                        support   = 1:10) {
  print("sim_model")
  print(sim_model)
  print("model")
  print(model)
  
  # 1. get simulation function
  sim <- get(sim_model)
  
  # 2. large sample for "truth"
  tmp0       <- sim(n = N0, setup = setup, cov = cov)
  data_large <- tmp0$dataframe
  cond_pmf   <- tmp0$cond_pmf    # only used if model != "np"
  cond_mean  <- tmp0$mean_time   # only used if model != "np"
  
  X_large <- data_large[, c("Delta","race","sex","cont")]
  
  # compute true conditional means & pmfs
  if(sim_model != "simulate_data_real"){
    true_mean_vec <- with(X_large,
                          mapply(cond_mean, race, sex, cont))
    true_pmf_mat  <- t(mapply(function(r,s,c) cond_pmf(r,s,c),
                              X_large$race,
                              X_large$sex,
                              X_large$cont))
  }else{ # for real simulation 
    true_mean_vec <- with(X_large,
                          mapply(cond_mean, Delta, race, sex, cont))
    true_pmf_mat  <- t(mapply(function(d,r,s,c) cond_pmf(d,r,s,c),
                              X_large$Delta,
                              X_large$race,
                              X_large$sex,
                              X_large$cont))
  }

  
  # true marginal pmf
  true_marg_pmf <- colMeans(true_pmf_mat)
  
  # storage
  se_vec      <- numeric(m)
  ise_vec     <- numeric(m)
  itv_vec     <- numeric(m)
  tv_marg_vec <- numeric(m)
  

  
  # 3. Monte Carlo replicates
  for(i in seq_len(m)) {
    # simulate smaller dataset & fit
    tmp    <- sim(n = n, setup = setup)
    df0    <- transform_df(tmp$dataframe, 11) 
    
    fit <- MIC(subset(df0, select = -T),
          user_formula = "cont + factor(sex) + factor(race)",
          method       = model,
          cov          = cov)
    
    # build estimated conditional pmf & mean on the large grid
    est_pmf_mat <- matrix(0, nrow = N0, ncol = length(support))
    est_mean_vec <- numeric(N0)
    
    for(j in seq_len(N0)) {
      if(model %in% c("roc", "pred", "pred.adj")){
        est_mean_vec[j] <- fit$meanT
      }else{
        pmfj <- pmf(fit$model,
                    cov      = list(race = X_large$race[j],
                                    sex   = X_large$sex[j],
                                    cont  = X_large$cont[j]),
                    max_time = 10)$y
        est_pmf_mat[j, ] <- pmfj
        est_mean_vec[j]  <- sum(support * pmfj)
      }
    }
    
    # SE 
    se_vec[i] <- (mean(true_mean_vec) - mean(est_mean_vec))**2 
    
    if(!(model %in% c("roc", "pred", "pred.adj"))){
      # ISE
      ise_vec[i] <- mean((true_mean_vec - est_mean_vec)^2)
      
      # ITV
      tv_cond    <- rowSums(abs(true_pmf_mat - est_pmf_mat)) / 2
      itv_vec[i] <- mean(tv_cond)
      
      # marginal TV
      est_marg_pmf   <- colMeans(est_pmf_mat)
      tv_marg_vec[i] <- sum(abs(est_marg_pmf - true_marg_pmf)) / 2
    }
  }
  
  if(model %in% c("roc", "pred", "pred.adj")){
    ise_vec <- NA
    itv_vec <- NA
    tv_marg_vec <- NA 
  }
  
  out <- list(
    model         = model,
    se            = se_vec, 
    ise           = ise_vec,      # ∫(E[T|X]-ĤE[T|X])² p(X)dX
    itv           = itv_vec,      # ∫TV(P(T|X),ĤP(T|X)) p(X)dX
    tv_marg       = tv_marg_vec,  # TV(E[ĤP(T|X)],P(T))
    true_marg_pmf = true_marg_pmf,
    n             = n,
    m             = m,
    setup         = setup
  )
  class(out) <- "mega"
  return(out)
}


#setup is fixed 
run_simulations <- function(setup, cov = NULL, m = 100, sim_models = c("simulate_data_po", "simulate_data_ph", "simulate_data_np", "simulate_data_real"), models = c("np", "ph", "po", "roc", "pred", "pred.adj")) {
  # Define the sample sizes and simulation model choices
  n_values <- c(100)
  sim_models <- sim_models
  models <- models
  
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
  
  # 1) Marginal mean squared error
  cat("MSE (marginal mean):", mean(x$se, na.rm = TRUE), "\n\n")
  
  # 2) Integrated squared error (conditional mean)
  if (all(is.na(x$ise))) {
    cat("ISE (conditional mean): Not available\n\n")
  } else {
    cat("ISE (conditional mean):", mean(x$ise, na.rm = TRUE), "\n\n")
  }
  
  # 3) Integrated TV distance (conditional)
  if (all(is.na(x$itv))) {
    cat("ITV (conditional TV): Not available\n\n")
  } else {
    cat("ITV (conditional TV):", mean(x$itv, na.rm = TRUE), "\n\n")
  }
  
  # 4) Marginal TV distance
  if (all(is.na(x$tv_marg))) {
    cat("Marginal TV distance: Not available\n")
  } else {
    cat("Marginal TV distance:", mean(x$tv_marg, na.rm = TRUE), "\n")
  }
  
  invisible(x)
}


print_outdated.mega <- function(x) {
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
  
  # Print the coverage error for the CI for the mean
  if(!all(is.na(x$coverage))) {
    # Compute the average coverage error: the proportion of replicates in which the CI did NOT cover the true mean.
    cov_error <- mean(x$coverage, na.rm = TRUE)
    cat("Coverage error (proportion of replicates where CI did NOT cover the true mean):", cov_error, "\n")
  } else {
    cat("Coverage error: Not available\n")
  }
  
  invisible(x)
}
