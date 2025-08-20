
library(parallel)

generic_sim <- function(n, m, setup, cov = NULL,
                        sim_model = "simulate_data_np",
                        model     = "np",
                        N0        = 10000,
                        support   = 1:10,
                        boot      = FALSE) {
  
  print("sim_model")
  print(sim_model)
  print("model")
  print(model)
  
  # 1. get simulation function
  sim <- get(sim_model)
  
  # 2. large sample for "truth"
  tmp0       <- sim(n = N0, setup = setup, cov = cov)
  data_large <- tmp0$dataframe
  cond_pmf   <- tmp0$cond_pmf
  cond_mean  <- tmp0$mean_time
  
  X_large <- data_large[, c("Delta","race","sex","cont")]
  
  if(sim_model != "simulate_data_real"){
    true_mean_vec <- with(X_large,
                          mapply(cond_mean, race, sex, cont))
    true_pmf_mat  <- t(mapply(function(r,s,c) cond_pmf(r,s,c),
                              X_large$race,
                              X_large$sex,
                              X_large$cont))
  } else {
    true_mean_vec <- with(X_large, mapply(cond_mean, Delta, race, sex, cont))
    true_pmf_mat  <- t(mapply(function(d,r,s,c) cond_pmf(d,r,s,c),
                              X_large$Delta,
                              X_large$race,
                              X_large$sex,
                              X_large$cont))
  }
  
  true_marg_pmf <- colMeans(true_pmf_mat)
  
  # 3. Monte Carlo replicates (parallel)
  results_list <- mclapply(seq_len(m), function(i) {
    tmp    <- sim(n = n, setup = setup)
    df0    <- transform_df(tmp$dataframe, 11) 
    
    fit <- MIC(subset(df0, select = -T),
               user_formula = "cont + factor(sex) + factor(race)",
               method       = model,
               cov          = cov, 
               boot         = boot)
    
    est_pmf_mat <- matrix(0, nrow = N0, ncol = length(support))
    est_mean_vec <- numeric(N0)
    
    for (j in seq_len(N0)) {
      if (model %in% c("roc", "pred", "pred.adj")) {
        est_mean_vec[j] <- fit$meanT
      } else {
        pmfj <- pmf(fit$model,
                    cov      = list(race = X_large$race[j],
                                    sex   = X_large$sex[j],
                                    cont  = X_large$cont[j]),
                    max_time = 10)$y
        est_pmf_mat[j, ] <- pmfj
        est_mean_vec[j]  <- sum(support * pmfj)
      }
    }
    
    true_mean <- mean(true_mean_vec)
    se_val <- (true_mean - mean(est_mean_vec))^2
    coverage_val <- as.numeric(true_mean >= fit$CI[1] & true_mean <= fit$CI[2])
    bias_val <- mean(est_mean_vec) - true_mean
    
    if (model %in% c("roc", "pred", "pred.adj")) {
      ise_val <- NA
      itv_val <- NA
      tv_marg_val <- NA
    } else {
      ise_val <- mean((true_mean_vec - est_mean_vec)^2)
      tv_cond <- rowSums(abs(true_pmf_mat - est_pmf_mat)) / 2
      itv_val <- mean(tv_cond)
      est_marg_pmf <- colMeans(est_pmf_mat)
      tv_marg_val  <- sum(abs(est_marg_pmf - true_marg_pmf)) / 2
    }
    
    list(se = se_val,
         coverage = coverage_val,
         bias = bias_val,
         ise = ise_val,
         itv = itv_val,
         tv_marg = tv_marg_val)
    
  }, mc.cores = 6) # number of cores 
  
  # Extract results from list
  se_vec      <- sapply(results_list, `[[`, "se")
  coverage    <- sapply(results_list, `[[`, "coverage")
  bias        <- sapply(results_list, `[[`, "bias")
  ise_vec     <- sapply(results_list, `[[`, "ise")
  itv_vec     <- sapply(results_list, `[[`, "itv")
  tv_marg_vec <- sapply(results_list, `[[`, "tv_marg")
  
  if (model %in% c("roc", "pred", "pred.adj")) {
    ise_vec <- NA
    itv_vec <- NA
    tv_marg_vec <- NA 
  }
  
  out <- list(
    model         = model,
    se            = se_vec, 
    bias          = bias,
    coverage      = coverage, 
    ise           = ise_vec,
    itv           = itv_vec,
    tv_marg       = tv_marg_vec,
    true_marg_pmf = true_marg_pmf,
    n             = n,
    m             = m,
    setup         = setup
  )
  class(out) <- "mega"
  return(out)
}


generic_sim_nonparell <- function(n, m, setup, cov = NULL,
                        sim_model = "simulate_data_np",
                        model     = "np",
                        N0        = 10000,
                        support   = 1:10,
                        boot      = FALSE) {
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
  #true_mean <- mean(data_large$T)
  
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
  bias        <- numeric(m)
  coverage    <- numeric(m)
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
          cov          = cov, 
          boot         = boot)
    
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
    
    # True mean
    true_mean <- mean(true_mean_vec)
    
    # SE 
    se_vec[i] <- (true_mean - mean(est_mean_vec))**2 
    
    # Coverage
    coverage[i] <- as.numeric(true_mean >= fit$CI[1] & true_mean <= fit$CI[2])
    
    # Bias
    bias[i] <- mean(est_mean_vec) - true_mean 
    
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
    bias          = bias,
    coverage      = coverage, 
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


# Setup is fixed 
run_simulations <- function(setup, cov = NULL, m = 100, boot = FALSE, sim_models = c("simulate_data_po", "simulate_data_ph", "simulate_data_np", "simulate_data_real"), models = c("np", "ph", "po", "roc", "pred", "pred.adj")) {
  # Define the sample sizes and simulation model choices
  n_values <- c(100, 500, 1000)
  sim_models <- sim_models
  models <- models
  
  # Store results for each combination
  results_list <- list()
  
  # Loop over each sample size and simulation/model choice
  for (n_val in n_values) {
    for (i in seq_along(sim_models)) {
      for (j in seq_along(models)) {
        sim_model <- sim_models[i]
        model <- models[j]
        
        # Call the generic simulation function
        sim_result <- generic_sim(n = n_val, m = m, setup = setup, cov = cov,
                                  sim_model = sim_model, model = model, boot = boot)
        
        # Create a key for the result list, e.g. "n_100_po"
        key <- paste("n", n_val, sim_model, "fit_model", model, sep = "_")
        results_list[[key]] <- sim_result 
      }
    }
  }
  
  return(results_list)
}


summarize_results <- function(results_list, file = "results.csv") {
  summary_rows <- list()
  
  for (key in names(results_list)) {
    res <- results_list[[key]]
    
    # parse name
    parts     <- strsplit(key, "_")[[1]]
    n_val     <- as.integer(parts[2])
    sim_model <- parts[3]
    fit_model <- parts[6]
    
    # warn if any NAs in the metric vectors
    for (metric in c("ise", "itv", "tv_marg")) {
      vec <- res[[metric]]
      if (!all(is.na(vec))) {
        prop_na <- mean(is.na(vec))
        if (prop_na > 0) {
          warning(sprintf(
            "%s: %.1f%% of %s replicates are NA",
            key, 100 * prop_na, toupper(metric)
          ), call. = FALSE)
        }
      }
    }
    
    # Summary
    mse      <- mean(res$se,      na.rm = TRUE)
    bias     <- mean(res$bias,    na.rm = TRUE)
    coverage <- mean(res$coverage, na.rm = TRUE)
    ise      <- if (all(is.na(res$ise)))     NA else mean(res$ise,     na.rm = TRUE)
    itv      <- if (all(is.na(res$itv)))     NA else mean(res$itv,     na.rm = TRUE)
    tv_marg  <- if (all(is.na(res$tv_marg))) NA else mean(res$tv_marg, na.rm = TRUE)
    
    this_row <- data.frame(
      n         = n_val,
      setup     = res$setup,
      sim_model = sim_model,
      fit_model = fit_model,
      MSE       = mse,
      Bias      = bias,
      Coverage  = coverage,
      ISE       = ise,
      ITV       = itv,
      TV_marg   = tv_marg,
      stringsAsFactors = FALSE
    )
    
    summary_rows[[key]] <- this_row
  }
  
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file = file, row.names = FALSE)
  message("Wrote summary to ", file)
  invisible(summary_df)
}


# Example usage:
# results_list <- run_simulations(setup = 0, m = 50)   
# summarize_results(results_list, file = "results.csv")

print.mega <- function(x) {
  cat(x$model, "Simulation Results\n")
  cat("------------------------\n")
  cat("Parameters: n =", x$n, ", m =", x$m, ", setup =", x$setup, "\n\n")
  
  # Marginal mean squared error
  cat("MSE (marginal mean):", mean(x$se, na.rm = TRUE), "\n\n")
  
  # Bias
  cat("Bias:", mean(x$bias, na.rm = TRUE), "\n\n")
  
  # Coverage
  cat("Coverage:", mean(x$coverage, na.rm = TRUE), "\n\n")
  
  # Integrated squared error (conditional mean)
  if (all(is.na(x$ise))) {
    cat("ISE (conditional mean): Not available\n\n")
  } else {
    cat("ISE (conditional mean):", mean(x$ise, na.rm = TRUE), "\n\n")
  }
  
  # Integrated TV distance (conditional)
  if (all(is.na(x$itv))) {
    cat("ITV (conditional TV): Not available\n\n")
  } else {
    cat("ITV (conditional TV):", mean(x$itv, na.rm = TRUE), "\n\n")
  }
  
  # Marginal TV distance
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
