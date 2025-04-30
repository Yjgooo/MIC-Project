simulate_data_real <- function(
    n,
    seed    = NULL,
    setup = 0,
    cov = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  if(setup == 1){
    # supports
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 0, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 5
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 0.8, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 0){ # makes T indep D | X; Will this make ph/po method win? 
    # supports
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 5
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 0, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 0.1){ # Makes T more strongly correlated with D even after conditioning on X; Will this make pred.adj win bigger?
    # supports
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 2
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 1, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
    
  }else if(setup == 0.2){
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 2
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 2, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 0.3){
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 2
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 3, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 0.4){
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 2
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 4, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 0.5){
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 2
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 1.5, sex = 0.4, race = -0.2, cont = -0.01)
    sd_T     = 2
  }else if(setup == 1.1){ # Makes sex to have strong effect on T
    # supports
    Delta_vals = -5:10
    T_vals     = 1:10
    # Delta ~ discrete normal conditional on X
    coefs_D = c(intercept = 1, sex = 0.5, race = 0, cont = 0.01)
    sd_D     = 5
    # T ~ discrete normal conditional on Delta and X
    coefs_T = c(intercept = 1, Delta = 0, sex = 3, race = -0.2, cont = -0.01)
    sd_T     = 2
  }
  
  # 1) Generate covariates
  sex  <- sample(0:1, n, TRUE)
  race <- sample(0:2, n, TRUE)
  cont <- rnorm(n, 65, 10)
  
  # 2) Simulate Delta
  Delta <- numeric(n)
  for (i in seq_len(n)) {
    mu_D <- coefs_D['intercept'] +
      coefs_D['sex']   * sex[i] +
      coefs_D['race']  * (race[i] - 1) +
      coefs_D['cont']  * cont[i]
    pmf_D <- dnorm(Delta_vals, mean = mu_D, sd = sd_D)
    pmf_D <- pmf_D / sum(pmf_D)
    Delta[i] <- sample(Delta_vals, 1, prob = pmf_D)
  }
  
  # 3) Simulate T
  T <- numeric(n)
  for (i in seq_len(n)) {
    mu_T <- coefs_T['intercept'] +
      coefs_T['Delta'] * Delta[i] +
      coefs_T['sex']   * sex[i]    +
      coefs_T['race']  * (race[i] - 1) +
      coefs_T['cont']  * cont[i]
    pmf_T <- dnorm(T_vals, mean = mu_T, sd = sd_T)
    pmf_T <- pmf_T / sum(pmf_T)
    T[i] <- sample(T_vals, 1, prob = pmf_T)
  }
  
  # 4) Indicator
  Indicator <- as.integer(Delta >= T)
  
  # 5) Assemble data.frame
  df <- data.frame(Delta = Delta, T = T, Indicator = Indicator,
                   sex = sex, race = race, cont = cont)
  
  # 6) Conditional PMF of T|X
  cond_pmf <- function(Delta, sex, race, cont) {
    mu_T <- coefs_T['intercept'] +
      coefs_T['Delta'] * Delta +
      coefs_T['sex']   * sex    +
      coefs_T['race']  * (race - 1) +
      coefs_T['cont']  * cont
    pmf <- dnorm(T_vals, mean = mu_T, sd = sd_T)
    pmf <- pmf / sum(pmf)
    names(pmf) <- paste0('t=', T_vals)
    pmf
  }
  
  # 7) Mean event time
  mean_time <- function(Delta, sex, race, cont) {
    pmf <- cond_pmf(Delta, sex, race, cont)
    sum(T_vals * pmf)
  }
  
  list(dataframe = df, cond_pmf = cond_pmf, mean_time = mean_time)
}


simulate_data_real_bad_news <- function(n, seed = NULL, setup = 0, cov = NULL,
                               coefs = c(intercept = 6, Delta = 0.2, sex = 0.5, race = -0.3)) {
  if (!is.null(seed)) set.seed(seed)
  
  # Support for Δ and T
  Delta_vals <- -10:10
  T_vals     <- 1:10
  
  # Pre-allocate
  Delta <- numeric(n)
  T     <- numeric(n)
  
  # Covariates
  sex  <- sample(0:1, n, replace = TRUE)
  race <- sample(0:2, n, replace = TRUE)
  age  <- rnorm(n, mean = 65, sd = 10)
  
  # Simulate Δ and T
  for (i in seq_len(n)) {
    # Baseline probabilities for Δ
    base_probs <- exp(-abs(Delta_vals) / 3)
    # Age adjustment
    age_adj <- if (age[i] > 65) exp(-abs(Delta_vals + 2) / 3) else rep(1, length(Delta_vals))
    Delta_probs <- base_probs * age_adj
    Delta_probs <- Delta_probs / sum(Delta_probs)
    Delta[i] <- sample(Delta_vals, 1, prob = Delta_probs)
    
    # Conditional mean for T using coefficients vector
    mean_T_i <- coefs["intercept"] +
      coefs["Delta"]     * Delta[i] +
      coefs["sex"]       * sex[i] +
      coefs["race"]      * (race[i] - 1)
    
    # Probabilities for T
    T_density <- dnorm(T_vals, mean = mean_T_i, sd = 1)
    T_probs   <- T_density / sum(T_density)
    T[i] <- sample(T_vals, 1, prob = T_probs)
  }
  
  # Indicator
  Indicator <- as.integer(Delta >= T)
  
  # Assemble data.frame
  df <- data.frame(
    Delta     = Delta,
    T         = T,
    Indicator = Indicator,
    sex       = sex,
    race      = race,
    cont      = age
  )
  
  # Conditional PMF of T | (Delta, sex, race, cont)
  cond_pmf <- function(Delta, sex, race, cont) {
    base_probs <- exp(-abs(Delta_vals) / 3)
    age_adj   <- if (cont > 65) exp(-abs(Delta_vals + 2) / 3) else rep(1, length(Delta_vals))
    Delta_probs <- base_probs * age_adj
    Delta_probs <- Delta_probs / sum(Delta_probs)
    
    mean_T <- coefs["intercept"] +
      coefs["Delta"] * Delta +
      coefs["sex"]   * sex +
      coefs["race"]  * (race - 1)
    
    T_density <- dnorm(T_vals, mean = mean_T, sd = 1)
    pmf <- T_density / sum(T_density)
    names(pmf) <- paste0("t=", T_vals)
    pmf
  }
  
  # Mean event time E[T | X]
  mean_time <- function(Delta, sex, race, cont) {
    pmf <- cond_pmf(Delta, sex, race, cont)
    sum(T_vals * pmf)
  }
  
  list(
    dataframe = df,
    cond_pmf   = cond_pmf,
    mean_time  = mean_time
  )
}



simulate_data_real_old <- function(n, seed = NULL, setup = 0, cov = NULL) { #setup is dummmy for now
  if (!is.null(seed)) set.seed(seed)
  
  # Define possible values:
  Delta_vals <- -10:10      # Δ corresponds to change in pain level
  T_vals     <- 1:10        # T corresponds to the patient's minimal important change
  
  # Pre-allocate vectors for Δ and T
  Delta <- numeric(n)
  T     <- numeric(n)
  
  # Covariate generation:
  sex  <- sample(0:1, n, replace = TRUE)         # 0/1 indicator for sex
  race <- sample(0:2, n, replace = TRUE)           # e.g., 0, 1, 2 for different race groups
  age  <- rnorm(n, mean = 65, sd = 10)             # age with mean 65 and sd 10
  
  # Simulate Δ (pain change) and T (minimal important change)
  for (i in 1:n) {
    # Create a baseline probability for Δ values: higher weight near 0 (smaller changes)
    base_probs <- exp(-abs(Delta_vals)/3)
    
    # Optionally, adjust Δ probabilities based on age: older patients might tend to show less improvement
    age_adj <- if (age[i] > 65) exp(-abs(Delta_vals + 2)/3) else rep(1, length(Delta_vals))
    Delta_probs <- base_probs * age_adj
    Delta_probs <- Delta_probs / sum(Delta_probs)
    
    # Sample Δ for the i-th patient:
    Delta[i] <- sample(Delta_vals, 1, prob = Delta_probs)
    
    # Now simulate T conditional on Δ and the covariates.
    # We set a mean for T that is influenced by Δ, sex and race.
    # For instance, a larger Δ (i.e. greater pain improvement) might be associated with a higher minimal change threshold.
    mean_T <- 6 + 0.2 * Delta[i] + 0.5 * sex[i] - 0.3 * (race[i] - 1)
    
    # Create probabilities for T using a normal kernel centered at mean_T.
    T_density <- dnorm(T_vals, mean = mean_T, sd = 1)
    T_probs <- T_density / sum(T_density)
    
    # Sample T for the i-th patient:
    T[i] <- sample(T_vals, 1, prob = T_probs)
  }
  
  # Compute the indicator: whether the patient perceives significant improvement
  # (i.e. the change in pain Δ meets or exceeds the minimal important change T)
  Indicator <- as.integer(Delta >= T)
  
  # Combine all variables into a data frame and return
  df <- data.frame(
    Delta     = Delta,
    T         = T,
    Indicator = Indicator,
    sex       = sex,
    race      = race,
    cont       = age
  )
  
  return(df)
}

run_simulations_real <- function(setup = 0, cov = NULL, m = 200) {
  # Define the sample sizes and simulation model choices
  n_values <- c(100, 500, 1000)
  sim_models <- c("simulate_data_real")
  models <- c("po", "ph", "np")
  
  # Container to store results for each combination
  results_list <- list()
  
  # Loop over each sample size and simulation/model choice
  for (n_val in n_values) {
    for (j in seq_along(models)) {
      sim_model <- sim_models[i]
      model <- ""
      
      # Call the generic simulation function
      sim_result <- generic_sim(n = n_val, m = m, setup = setup, cov = cov,
                                sim_model = sim_model, model = model)
      
      # Create a key for the result list, e.g. "n_100_po"
      key <- paste("n", n_val, sim_model, "fit_model", model, sep = "_")
      results_list[[key]] <- sim_result 
    }
  }
  
  return(results_list)
}
