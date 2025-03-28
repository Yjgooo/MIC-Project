simulate_data_real <- function(n, seed = NULL, setup = 0, cov = NULL) { #setup is dummmy for now
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
    # For instance, a larger Δ (i.e. greater pain improvement) might be associated with a lower minimal change threshold.
    mean_T <- 6 - 0.2 * Delta[i] + 0.5 * sex[i] - 0.3 * (race[i] - 1)
    
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
