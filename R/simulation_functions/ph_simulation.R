library(fastDummies)
library(simsurv)

simulate_data_ph <- function(n, seed = NULL, setup = 0, cov = NULL) {
  #— 0. Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  #— 1. Simulate covariates
  race <- factor(
    sample(0:2, size = n, replace = TRUE, prob = c(0.4, 0.2, 0.4)),
    levels = c(0, 1, 2)
  )
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
  
  #— 3. Baseline hazards and regression coefficients
  baseline_hazard <- c(
    0.01, 0.03, 0.10, 0.20, 0.30,
    0.35, 0.40, 0.50, 0.60, 0.04
  )
  betas <- switch(as.character(setup),
                  "0" = c(
                    race_1 = -0.5,
                    race_2 =  1.0,
                    sex_1  =  0.2,
                    cont   =  0.01
                  ),
                  stop("Invalid setup")
  )
  
  #— 4. Override covariates if user-supplied
  if (!is.null(cov)) {
    for (nm in names(cov)) {
      if (nm %in% names(cov_data_sim)) {
        cov_data_sim[[nm]] <- rep(cov[[nm]], n)
      } else {
        stop("Covariate ", nm, " not found.")
      }
    }
  }
  
  #— 5. Compute linear predictor
  lp <- as.numeric(as.matrix(cov_data_sim[, names(betas)]) %*% betas)
  
  #— 6. Sampler for discrete‐time PH model
  sample_discrete_ph <- function(lp_i) {
    for (t in seq_along(baseline_hazard)) {
      p_t <- pmin(baseline_hazard[t] * exp(lp_i), 1)
      if (rbinom(1, 1, p_t) == 1) return(t)
    }
    length(baseline_hazard)
  }
  
  #— 7. Draw event times T and Delta
  T_val  <- sapply(lp, sample_discrete_ph)
  Delta  <- sapply(lp, sample_discrete_ph)
  
  #— 8. Build final data frame
  final_data <- merge(
    cov_data,
    cov_data_sim[, setdiff(names(cov_data_sim), "cont")],
    by = "id"
  )
  final_data$T         <- T_val
  final_data$Delta     <- Delta
  final_data$Indicator <- as.integer(Delta >= T_val)
  
  #— 9. Corrected conditional PMF: P(T = t | X)
  cond_pmf <- function(race, sex, cont) {
    # coerce factors to integers
    race <- as.integer(as.character(race))
    sex  <- as.integer(as.character(sex))
    if (!race %in% 0:2) stop("race must be 0,1,2")
    if (!sex  %in% 0:1) stop("sex must be 0,1")
    
    # rebuild linear predictor
    lp_i <- betas["race_1"] * as.numeric(race == 1) +
      betas["race_2"] * as.numeric(race == 2) +
      betas["sex_1"]  * as.numeric(sex  == 1) +
      betas["cont"]   * cont
    
    # hazards and survival
    n_t   <- length(baseline_hazard)
    p_vec <- pmin(baseline_hazard * exp(lp_i), 1)
    q_vec <- 1 - p_vec
    
    # build pmf with all remaining mass at t = n_t
    pmf    <- numeric(n_t)
    S_prev <- 1
    for (t in seq_len(n_t - 1)) {
      pmf[t]  <- S_prev * p_vec[t]
      S_prev  <- S_prev * q_vec[t]
    }
    pmf[n_t] <- S_prev           # collapse tail mass here
    
    names(pmf) <- paste0("t=", seq_len(n_t))
    pmf
  }
  
  #— 10. Mean event time from the PMF
  mean_time <- function(race, sex, cont) {
    pmf <- cond_pmf(race, sex, cont)
    sum(seq_along(pmf) * pmf)
  }
  
  #— 11. Return
  list(
    dataframe = final_data,
    cond_pmf  = cond_pmf,
    mean_time = mean_time
  )
}


simulate_data_ph_wrong_cond_pm <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # simulate covariates
  race <- factor(
    sample(0:2, size = n, replace = TRUE, prob = c(0.4, 0.2, 0.4)),
    levels = c(0, 1, 2)
  )
  sex  <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  cov_data <- data.frame(id = 1:n, race = race, sex = sex, cont = cont)
  
  # one‑hot encode
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns         = c("race", "sex"),
    remove_first_dummy     = TRUE,
    remove_selected_columns = TRUE
  )
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  
  # baseline hazards and betas (no intercept)
  baseline_hazard <- c(0.01, 0.03, 0.10, 0.20,
                       0.30, 0.35, 0.40, 0.50,
                       0.60, 0.04)
  
  betas <- switch(as.character(setup),
                  "0" = c(
                    race_1 = -0.5,
                    race_2 =  1.0,
                    sex_1  =  0.2,
                    cont   =  0.01
                  ),
                  stop("Invalid setup")
  )
  
  # override covariates if requested
  if (!is.null(cov)) {
    for (nm in names(cov)) {
      if (nm %in% names(cov_data_sim)) {
        cov_data_sim[[nm]] <- rep(cov[[nm]], n)
      } else {
        stop("Covariate ", nm, " not found.")
      }
    }
  }
  
  # linear predictor
  lp <- as.numeric(as.matrix(cov_data_sim[, names(betas)]) %*% betas)
  
  # sampling function for PH (truncated support 1:10)
  sample_discrete_ph <- function(lp_i) {
    for (t in seq_along(baseline_hazard)) {
      p_t <- pmin(baseline_hazard[t] * exp(lp_i), 1)
      if (rbinom(1, 1, p_t) == 1) return(t)
    }
    length(baseline_hazard)
  }
  
  # draw event times
  T_val <- sapply(lp, sample_discrete_ph)
  Delta <- sapply(lp, sample_discrete_ph)
  
  # assemble data
  final_data <- merge(
    cov_data,
    cov_data_sim[, setdiff(names(cov_data_sim), "cont")],
    by = "id"
  )
  final_data$T         <- T_val
  final_data$Delta     <- Delta
  final_data$Indicator <- as.integer(Delta >= T_val)
  
  # conditional PMF: P(T = t | X)
  cond_pmf <- function(race, sex, cont) {
    if (!race %in% 0:2) stop("race must be 0,1,2")
    if (!sex  %in% 0:1) stop("sex must be 0,1")
    lp_i <- betas["race_1"] * as.numeric(race == 1) +
      betas["race_2"] * as.numeric(race == 2) +
      betas["sex_1"]  * as.numeric(sex == 1) +
      betas["cont"]   * cont
    p_vec <- pmin(baseline_hazard * exp(lp_i), 1)
    q_vec <- 1 - p_vec
    S_prev <- 1
    pmf <- numeric(length(baseline_hazard))
    for (t in seq_along(baseline_hazard)) {
      pmf[t] <- S_prev * p_vec[t]
      S_prev <- S_prev * q_vec[t]
    }
    names(pmf) <- paste0("t=", seq_along(pmf))
    pmf
  }
  
  # mean event time from the PMF
  mean_time <- function(race, sex, cont) {
    pmf <- cond_pmf(race, sex, cont)
    sum(seq_along(pmf) * pmf)
  }
  
  list(
    dataframe = final_data,
    cond_pmf  = cond_pmf,
    mean_time = mean_time
  )
}



ph_sim <- function(n,m, setup, cov = NULL){ 
  
  #this can be either 
  data_large <- simulate_data_ph(n = 10000, setup = setup, cov =cov)
  true_mean <- mean(data_large$T)
  
  full_support <- c(1:10)
  true_pmf <- prop.table(table(factor(data_large$T, levels = full_support)))
  
  
  vec <- numeric(m)
  tv_dist_vec = numeric(m)
  
  # models on restricted data;
  for(i in 1:m){
    data_oracle <- simulate_data_ph(n = n, setup = setup)
    data_observed <- subset(data_oracle, select = -T) #dataframe without T
    df_ <- transform_df(data_observed, score_max = 11)
    
    fit_ph <- MIC(df_, user_formula = "cont + factor(sex) + factor(race)","ph", cov)
    print("*********")
    print(true_pmf)
    print(pmf(fit_ph$model, max_time = 10)$y)
    
    vec[i] <- abs(true_mean - fit_ph$meanT)
    tv_dist_vec[i] <- tv_distance(true_pmf, pmf(fit_ph$model, cov, max_time = 10)$y)
  }
  
  result <- list(
    abs_dif = vec,
    tv_dist = tv_dist_vec, #abs_dif and tv_dist are vectors 
    n = n,
    m = m,
    setup = setup
  )
  
  class(result) <- "ph_sim"
  return(result) #will need to skips the nulls (caused by Max iteration reached)
}

print.ph_sim <- function(x) {
  cat("ph Simulation Results\n")
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