library(fastDummies)
library(simsurv)

simulate_data_ph <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  race <- factor(sample(0:2, size = n, replace = TRUE, prob = c(0.4, 0.2, 0.4)), levels = c(0, 1, 2))
  sex <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  
  # Create raw covariate dataframe
  cov_data <- data.frame(
    id = 1:n,
    race = race,
    sex = sex,
    cont = cont
  )
  
  # One-hot encode race and sex
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns = c("race", "sex"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
  
  # Combine with continuous variable
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  cov_data_sim$intercept <- 1
  
  # Strong effect of race for large mean survival difference
  baseline_hazard <- c(0.01, 0.03, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.04)
  
  switch(as.character(setup),
         "0" = { 
           betas <- c(intercept = 1, race_1 = -0.5, race_2 = 1, sex_1 = 0.2, cont = 0.01)
           #betas <- c(intercept = 1.2, race_1 = -1.5, race_2 = 0.5, sex_1 = 0.2, cont = 0.01)
         },
         stop("Invalid setup value")
  )
  
  # Override covariate values if given
  if (!is.null(cov)) {
    for (name in names(cov)) {
      if (name %in% names(cov_data_sim)) {
        cov_data_sim[[name]] <- rep(cov[[name]], n)
      } else {
        stop(paste("Covariate", name, "not found in the data."))
      }
    }
  }
  
  # Linear predictor for PH model
  linpred <- as.matrix(cov_data_sim[, names(betas)]) %*% betas
  
  # Sampling function for PH model (truncated to ensure support {1,...,10})
  sample_discrete_ph <- function(lp, baseline_hazard) {
    for (t in seq_along(baseline_hazard)) {
      hazard_t <- baseline_hazard[t] * exp(lp)
      hazard_t <- pmin(hazard_t, 1)
      event <- rbinom(1, 1, hazard_t)
      if (event == 1) return(t)
    }
    return(length(baseline_hazard))  # Truncate to 10 instead of 11
  }
  
  # Sample T and Delta independently with the same conditional distribution
  T_val <- sapply(linpred, sample_discrete_ph, baseline_hazard = baseline_hazard)
  Delta <- sapply(linpred, sample_discrete_ph, baseline_hazard = baseline_hazard)
  
  # Final output
  final_data <- merge(cov_data, cov_data_sim[, !(names(cov_data_sim) %in% c("cont"))], by = "id")
  final_data$T <- T_val
  final_data$Delta <- Delta
  final_data$Indicator <- as.integer(final_data$Delta >= final_data$T)
  
  return(final_data)
}



#Tried to make the correlation between T and Delta strong
simulate_data_ph_1 <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  race <- factor(sample(0:2, size = n, replace = TRUE), levels = c(0, 1, 2))
  sex <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  
  # Create a data frame with the raw covariates and an identifier
  cov_data <- data.frame(
    id = 1:n,
    race = race,
    sex = sex,
    cont = cont
  )
  
  # One-hot encode the categorical variables using fastDummies
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns = c("race", "sex"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
  
  # Combine the one-hot encoded variables with id and continuous variable
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  
  # Specify the baseline hazard for the event time T (times 1 to 10)
  baseline_hazard <- c(0.01, 0.03, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.04)
  
  # Set coefficients for the linear predictors for T and for Delta
  switch(as.character(setup),
         "0" = { 
           betas <- c(race_1 = 0.5, race_2 = -0.3, sex_1 = 0.2, cont = 0.01)
           gammas <- c(race_1 = 0.3, race_2 = -0.2, sex_1 = 0.1, cont = 0.005)
         },
         stop("Invalid setup value")
  )
  
  # If fixed covariate values are provided, override the generated values
  if (!is.null(cov)) {
    for (name in names(cov)) {
      if (name %in% names(cov_data_sim)) {
        cov_data_sim[[name]] <- rep(cov[[name]], n)
      } else {
        stop(paste("Covariate", name, "not found in the data."))
      }
    }
  }
  
  # Compute the linear predictor for T and generate event times using a proportional hazards model
  linpred <- as.matrix(cov_data_sim[, names(betas)]) %*% betas
  
  # Modified function: use exp(-lp) so that higher lp (and thus higher covariate values)
  # lead to lower hazards and therefore later (larger) event times T.
  sample_discrete_ph <- function(lp, baseline_hazard) {
    for (t in seq_along(baseline_hazard)) {
      hazard_t <- baseline_hazard[t] * exp(-lp)
      hazard_t <- pmin(hazard_t, 1)
      event <- rbinom(1, 1, hazard_t)
      if (event == 1) return(t)
    }
    return(length(baseline_hazard) + 1)
  }
  
  event_times <- sapply(linpred, sample_discrete_ph, baseline_hazard)
  
  # Generate Delta as a discrete outcome.
  # Change the support from -10:10 to 0:10 so that Delta is always positive.
  possible_delta <- 0:10
  sigma <- 2  # Spread parameter (smaller sigma concentrates probability more strongly at 3)
  # Baseline log-probabilities: highest at 3 and falling off symmetrically
  eta_delta <- -((possible_delta - 3)^2) / (2 * sigma^2)
  
  # Compute the linear predictor for Delta from the covariates using gammas
  lp_delta <- as.vector(as.matrix(cov_data_sim[, names(gammas)]) %*% gammas)
  
  Delta <- sapply(1:n, function(i) {
    # Incorporate the covariate effect by shifting the baseline probabilities
    unnormalized <- exp(eta_delta + possible_delta * lp_delta[i])
    probs <- unnormalized / sum(unnormalized)
    sample(possible_delta, size = 1, prob = probs)
  })
  
  # Merge the simulation output with the original covariate data
  final_data <- merge(cov_data, cov_data_sim[, !(names(cov_data_sim) %in% c("cont"))], by = "id")
  final_data$T <- pmin(event_times, 10)
  final_data$Delta <- Delta
  final_data$Indicator <- as.integer(final_data$Delta >= final_data$T)
  
  return(final_data)
}


simulate_data_ph_old <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  race <- factor(sample(0:2, size = n, replace = TRUE), levels = c(0, 1, 2))
  sex <- factor(sample(0:1, size = n, replace = TRUE), levels = c(0, 1))
  cont <- rnorm(n, mean = 80, sd = sqrt(5))
  
  eta <- c(-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 3, 2.5, 2.5, 2.5, 2, 1, 1, 0.5, 0.5, 0.5, 0.3)
  D_probs <- exp(eta) / sum(exp(eta))
  Delta <- sample(-10:10, size = n, replace = TRUE, prob = D_probs)
  
  cov_data <- data.frame(
    id = 1:n,
    race = race,
    sex = sex,
    cont = cont,
    Delta = Delta
  )
  
  cov_data_dummy <- fastDummies::dummy_cols(
    cov_data,
    select_columns = c("race", "sex"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
  
  cov_data_sim <- cov_data_dummy[, c("id", "race_1", "race_2", "sex_1", "cont")]
  
  # Specify the baseline discrete hazard (hand-specified for times 1 to 10)
  baseline_hazard <- c(0.01, 0.03, 0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.04) #latter values needs to be larger...
  
  switch(as.character(setup),
         "0" = { 
           betas <- c(race_1 = 0.5, race_2 = -0.3, sex_1 = 0.2, cont = 0.01)
         },
         stop("Invalid setup value")
  )
  
  if(!is.null(cov)){
    for (name in names(cov)) {
      if (name %in% names(cov_data_sim)) {
        cov_data_sim[[name]] <- rep(cov[[name]], n)
      } else {
        stop(paste("Covariate", name, "not found in the data."))
      }
    }
  }
  
  linpred <- as.matrix(cov_data_sim[, names(betas)]) %*% betas
  
  sample_discrete_ph <- function(lp, baseline_hazard) {
    for (t in seq_along(baseline_hazard)) {
      hazard_t <- baseline_hazard[t] * exp(lp)
      hazard_t <- pmin(hazard_t, 1)
      event <- rbinom(1, 1, hazard_t)
      if (event == 1) return(t)
    }
    return(length(baseline_hazard) + 1)
  }
  
  event_times <- sapply(linpred, sample_discrete_ph, baseline_hazard)
  
  final_data <- merge(cov_data, cov_data_sim[, !(names(cov_data_sim) %in% c("cont"))], by = "id")
  final_data$T <- pmin(event_times, 10)
  final_data$Indicator <- as.integer(final_data$Delta >= final_data$T)
  
  return(final_data)
}


ph_sim <- function(n,m, setup, cov = NULL){ 
  
  #this can be either 
  data_large <- simulate_data_ph(n = 10000, setup = setup, cov =cov)
  true_mean <- mean(data_large$T)
  
  full_support <- c(1:10)
  true_pmf <- prop.table(table(factor(data_large$T, levels = full_support)))
  
  
  vec <- numeric(m)
  tv_dist_vec = numeric(m)
  
  #models on restricted data;
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