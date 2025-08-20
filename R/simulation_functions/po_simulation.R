library(fastDummies)

simulate_data_po <- function(n, seed = NULL, setup = 0, cov = NULL) {
  #— set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  #— 1. Generate covariates
  race <- factor(sample(0:2, size = n, replace = TRUE), levels = c(0, 1, 2))
  sex  <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  
  cov_data <- data.frame(id = 1:n, race = race, sex = sex, cont = cont)
  
  #— 2. One‑hot encode race and sex (drop first level)
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns         = c("race", "sex"),
    remove_first_dummy     = TRUE,
    remove_selected_columns = TRUE
  )
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  
  #— 3. Baseline odds and coefficients
  baseline_odds <- c(
    0.01, 0.03, 0.10, 0.20, 0.30,
    0.35, 0.40, 0.50, 0.60, 0.04
  )
  betas <- switch(as.character(setup),
                  "0" = c(race_1 =  2.5,
                          race_2 = -2.5,
                          sex_1  =  0.2,
                          cont   =  0.01),
                  stop("Invalid setup value")
  )
  
  #— 4. Override covariates if user-supplied values are given
  if (!is.null(cov)) {
    for (nm in names(cov)) {
      if (nm %in% names(cov_data_sim)) {
        cov_data_sim[[nm]] <- rep(cov[[nm]], n)
      } else {
        stop("Covariate ", nm, " not found.")
      }
    }
  }
  
  #— 5. Compute linear predictor for each subject
  linpred <- as.numeric(as.matrix(cov_data_sim[, names(betas)]) %*% betas)
  
  #— 6. Sampler for the discrete proportional‐odds model
  sample_discrete_po <- function(lp_i) {
    for (t in seq_along(baseline_odds)) {
      odds_t <- baseline_odds[t] * exp(lp_i)
      p_t    <- odds_t / (1 + odds_t)
      if (rbinom(1, 1, p_t) == 1) return(t)
    }
    length(baseline_odds)
  }
  
  #— 7. Draw event times T and Delta
  T_val  <- sapply(linpred, sample_discrete_po)
  Delta  <- sapply(linpred, sample_discrete_po)
  
  #— 8. Merge back to original covariate frame
  final_data <- merge(
    cov_data,
    cov_data_sim[, setdiff(names(cov_data_sim), "cont")],
    by = "id"
  )
  final_data$T         <- pmin(T_val, length(baseline_odds))
  final_data$Delta     <- Delta
  final_data$Indicator <- as.integer(Delta >= T_val)
  
  #— 9. Conditional PMF function P(T = t | X)
  cond_pmf <- function(race, sex, cont) {
    # coerce factors to integers if needed
    race <- as.integer(as.character(race))
    sex  <- as.integer(as.character(sex))
    if (!race %in% 0:2) stop("race must be 0,1,2")
    if (!sex  %in% 0:1) stop("sex must be 0,1")
    
    # rebuild linear predictor
    lp_i <- betas["race_1"] * as.numeric(race == 1) +
      betas["race_2"] * as.numeric(race == 2) +
      betas["sex_1"]  * as.numeric(sex  == 1) +
      betas["cont"]   * cont
    
    n_t   <- length(baseline_odds)
    odds  <- baseline_odds * exp(lp_i)
    p_vec <- odds / (1 + odds)
    q_vec <- 1 - p_vec
    
    # build pmf so that all remaining mass at t = n_t
    pmf    <- numeric(n_t)
    S_prev <- 1
    for (t in seq_len(n_t - 1)) {
      pmf[t]  <- S_prev * p_vec[t]
      S_prev  <- S_prev * q_vec[t]
    }
    pmf[n_t] <- S_prev
    
    names(pmf) <- paste0("t=", seq_len(n_t))
    pmf
  }
  
  #— 10. Mean event time from the pmf
  mean_time <- function(race, sex, cont) {
    pmf <- cond_pmf(race, sex, cont)
    sum(seq_along(pmf) * pmf)
  }
  
  #— Return the data frame plus helper functions
  list(
    dataframe = final_data,
    cond_pmf   = cond_pmf,
    mean_time  = mean_time
  )
}



simulate_data_po_old <- function(n, seed = NULL, setup = 0, cov = NULL) {
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