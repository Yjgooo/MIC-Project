library(fastDummies)

simulate_data_po <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  race <- factor(sample(0:2, size = n, replace = TRUE), levels = c(0, 1, 2))
  sex  <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  
  # Create a data frame with the raw covariates and an identifier
  cov_data <- data.frame(
    id   = 1:n,
    race = race,
    sex  = sex,
    cont = cont
  )
  
  # One-hot encode categorical variables using fastDummies
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns = c("race", "sex"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
  
  # Combine one-hot encoded variables with id and continuous variable
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  
  # Baseline odds for the event time T (times 1 to 10)
  baseline_odds_T <- c(0.01, 0.03, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.04)
  
  # Set coefficients for T (and implicitly Delta) linear predictor
  switch(as.character(setup),
         "0" = { 
           betas <- c(race_1 = 2.5,  race_2 = -2.5, sex_1 = 0.2, cont = 0.01)
         },
         stop("Invalid setup value")
  )
  
  # Override covariate values if fixed covariates are provided
  if (!is.null(cov)) {
    for (name in names(cov)) {
      if (name %in% names(cov_data_sim)) {
        cov_data_sim[[name]] <- rep(cov[[name]], n)
      } else {
        stop(paste("Covariate", name, "not found in the data."))
      }
    }
  }
  
  # --- Linear predictor ---
  linpred <- as.matrix(cov_data_sim[, names(betas)]) %*% betas
  
  sample_discrete_po <- function(lp, baseline_odds) {
    for (t in seq_along(baseline_odds)) {
      odds_t <- baseline_odds[t] * exp(lp)
      p_event <- odds_t / (1 + odds_t)
      p_event <- pmin(p_event, 1)
      event <- rbinom(1, 1, p_event)
      if (event == 1) return(t)
    }
    # Instead of returning 11, just return 10
    return(length(baseline_odds))
  }
  
  # --- Sample T ---
  event_times <- sapply(linpred, sample_discrete_po, baseline_odds = baseline_odds_T)
  T_val <- pmin(event_times, 10)
  
  # --- Sample Delta using the same conditional distribution as T ---
  Delta <- sapply(linpred, sample_discrete_po, baseline_odds = baseline_odds_T)
  
  # --- Final data assembly ---
  final_data <- merge(cov_data, cov_data_sim[, !(names(cov_data_sim) %in% c("cont"))], by = "id")
  final_data$T <- T_val
  final_data$Delta <- Delta
  final_data$Indicator <- as.integer(final_data$Delta >= final_data$T)
  
  return(final_data)
}

po_sim <- function(n, m, setup, cov = NULL){ 
  # Use a large sample to get the true distribution of T
  data_large <- simulate_data_po(n = 10000, setup = setup, cov = cov)
  true_mean <- mean(data_large$T)
  
  full_support <- c(1:10)
  true_pmf <- prop.table(table(factor(data_large$T, levels = full_support)))
  
  vec <- numeric(m)
  tv_dist_vec <- numeric(m)
  
  # Models on restricted data
  for(i in 1:m){
    data_oracle <- simulate_data_po(n = n, setup = setup)
    data_observed <- subset(data_oracle, select = -T) # dataframe without T
    df_ <- transform_df(data_observed, score_max = 11)
    
    fit_po <- MIC(df_, user_formula = "cont + factor(sex) + factor(race)", "po", cov)
    
    vec[i] <- abs(true_mean - fit_po$meanT)
    tv_dist_vec[i] <- tv_distance(true_pmf, pmf(fit_po$model, cov, max_time = 10)$y)
  }
  
  result <- list(
    abs_dif = vec,
    tv_dist = tv_dist_vec, # abs_dif and tv_dist are vectors 
    n = n,
    m = m,
    setup = setup
  )
  
  class(result) <- "po_sim"
  return(result) # will need to skips the nulls (caused by Max iteration reached)
}

print.po_sim <- function(x) {
  cat("po Simulation Results\n")
  cat("------------------------\n")
  cat("Parameters: n =", x$n, ", m =", x$m, ", setup =", x$setup, "\n\n")
  cat("MSE:", mean((x$abs_dif)**2, na.rm = TRUE), "\n")
  cat("95th percentile of MSE:", 
      as.numeric(quantile((x$abs_dif)**2, probs = 0.95, na.rm = TRUE)), "\n")
  cat("MAD:", mean(x$abs_dif, na.rm = TRUE), "\n")
  cat("95th percentile of MADs:", 
      as.numeric(quantile(x$abs_dif, probs = 0.95, na.rm = TRUE)), "\n")
  cat("Average TV distance:", mean(x$tv_dist, na.rm = TRUE), "\n")
  
  invisible(x)
}