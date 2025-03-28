
simulate_data_np <- function(n, seed = NULL, setup = 0, cov = NULL) { #cov = NULL is for facilitating generic simulation
  
  k <- 10 #so score goes from 1 to 5
  
  # Define possible values for Delta and T
  delta_vals <- -k:k   # Delta takes values from -4 to 4 (9 categories)
  T_vals     <- 1:k    # T takes values from 1 to 5 (5 categories)
  
  # Pre-allocate vectors for Delta and T.
  max_score <- numeric(1)
  Delta <- numeric(n)
  T <- numeric(n)
  trueT <- numeric(1)
  
  #Covariates generation 
  sex <- sample(0:1, n, replace = TRUE, prob = rep(1/2, 2))
  race <- sample(0:2, n, replace = TRUE, prob = rep(1/3, 3))
  age <- rnorm(n, mean = 80, sd = 5)
  
  #Sample Delta and T conditionally on the covariates
  # Sample Delta and T conditionally on the covariates using switch
  switch(as.character(setup),
         "0" = { #Uniform T and delta 
           
           for (i in 1:n) {
             weights_delta <- rep(1/(2*k + 1), 2*k + 1)
             Delta[i] <- sample(delta_vals, 1, prob = weights_delta)
             weights_T <- rep(1/k, k)
             T[i] <- sample(T_vals, 1, prob = weights_T)
           }
         },
         "1" = { #(X,D,T) mutually independent; allow for specification of marginal distribution of T and D 
           
           T_probs <- exp(c(2, 1, 0, 0))
           T_probs <- T_probs / sum(T_probs)
           
           eta_ <- c(0, 0.1, 0.2, 0.5, 3, 2, 1, 0.5, 0.1)
           D_probs <- exp(eta_) / sum(exp(eta_))
           
           T <- sample(T_vals, size = n, prob = T_probs, replace = TRUE)
           Delta <- sample(delta_vals, n, prob = D_probs, replace = TRUE)
           
         },
         "2" = { #T and D are not independent; T is inversely related to D; violation of the independence assumption 
           
           d_vec <- c(0, 0.1, 0.2, 0.5, 3, 2, 1, 0.5, 0.1)
           D_probs <- exp(d_vec) / sum(exp(d_vec))
           Delta <- sample(delta_vals, n, prob = D_probs, replace = TRUE)
           
           T <- pmin(4 - Delta + 1, 4) #max_score = 4
         },
         
         "3" = { 
           
           eta_ <- c(0, 0.1, 0.2, 0.5, 3, 2, 1, 0.5, 0.1)
           
           vec1 <- c(2, 1, 0, 0)
           vec2 <- c(0, 1, 2, 0)
           vec3 <- c(0, 3, 1, 0)
           
           mat_race <- cbind(vec1, vec2, vec3)
           
           vec1_ <- c(0.3, 1, -0.3, -0.1)
           vec2_ <- c(0, 0, 1, 0.3)
           
           mat_sex <- cbind(vec1_, vec2_)
           
           for (i in 1:n) {
             T_probs <- exp(mat_race[, race[i] + 1] + mat_sex[, sex[i] + 1])
             T_probs <- T_probs / sum(T_probs)
             
             D_probs <- exp(eta_)
             D_probs <- D_probs / sum(D_probs)
             
             T[i] <- sample(T_vals, size = 1, prob = T_probs)
             Delta[i] <- sample(delta_vals, 1, prob = D_probs)
           }
           
         },
         
         stop("Invalid setup value")
  )
  
  #Compute the indicator variable (feeling significantly better)
  Indicator <- as.integer(Delta >= T) #indicator whether patient feel significantly improved 
  
  #Combine into a data frame
  # 'race' remains as a numeric variable.
  df <- data.frame(
    Delta = Delta,
    T = T,
    Indicator = Indicator,
    race = race,
    sex = sex,
    cont = age
  )
  return(df)
}


np_sim <- function(n,m, setup, cov = NULL){ #cov = NULL is a just a placeholder to faciliate generic simulation; never used in the function
  
  #truth simulation 
  vec <- numeric(m)
  data_large <- simulate_data_np(10000, setup)
  true_mean <- mean(data_large$T)
  
  
  full_support <- c(1:10)
  true_pmf <- prop.table(table(factor(data_large$T, levels = full_support)))
  
  #repeated experiment
  tv_dist_vec = numeric(m)
  
  for(i in 1:m){
    data_oracle <- simulate_data_np(n, setup)
    data_observed <- subset(data_oracle, select = -T) #dataframe without T
    df <- transform_df(data_observed, score_max = 11)
    
    true_mean = mean(data_oracle$T)
    fit_np <- MIC(df, "np")
    vec[i] <- abs(true_mean - fit_np$meanT)
    #print("**********")
    #print(true_pmf)
    #print(pmf(fit_np$model, max_time = 10)$y)
    tv_dist_vec[i] <- tv_distance(true_pmf, pmf(fit_np$model, max_time = 10)$y) #NPMLE simulation only considers max_time = 4 
  }
  
  result <- list(
    abs_dif = vec,
    tv_dist = mean(tv_dist_vec),
    n = n,
    m = m,
    setup = setup
  )
  
  class(result) <- "np_sim"
  return(result)
}

print.np_sim <- function(x, ...) {
  cat("NPMLE Simulation Results\n")
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