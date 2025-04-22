
simulate_data_np <- function(n, seed = NULL, setup = 0, cov = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  k <- 10
  delta_vals <- 1:k
  T_vals     <- 1:k
  
  # simulate covariates
  sex  <- sample(0:1, n, replace = TRUE)
  race <- sample(0:2, n, replace = TRUE)
  age  <- rnorm(n, mean = 80, sd = 5)
  
  # draw Delta and T
  Delta <- integer(n)
  T     <- integer(n)
  
  switch(as.character(setup),
         "0" = {
           for (i in seq_len(n)) {
             Delta[i] <- sample(delta_vals, 1, prob = rep(1/k, k))
             T[i]     <- sample(T_vals, 1,     prob = rep(1/k, k))
           }
         },
         "1" = {
           T_probs <- exp(c(2,1,0,0)); T_probs <- T_probs/sum(T_probs)
           Delta_probs <- exp(c(0,0.1,0.2,0.5,3,2,1,0.5,0.1))
           Delta_probs <- Delta_probs/sum(Delta_probs)
           T     <- sample(T_vals, n, prob = rep_len(T_probs, k),     replace = TRUE)
           Delta <- sample(delta_vals, n, prob = Delta_probs,         replace = TRUE)
         },
         "2" = {
           Delta_probs <- exp(c(0,0.1,0.2,0.5,3,2,1,0.5,0.1))
           Delta_probs <- Delta_probs/sum(Delta_probs)
           Delta <- sample(delta_vals, n, prob = Delta_probs, replace = TRUE)
           T     <- pmin(4 - Delta + 1, 4)
         },
         "3" = {
           eta_ <- c(0,0.1,0.2,0.5,3,2,1,0.5,0.1)
           vec1  <- c(2,1,0,0); vec2 <- c(0,1,2,0); vec3 <- c(0,3,1,0)
           mat_race <- cbind(vec1, vec2, vec3)
           vec1_ <- c(0.3,1,-0.3,-0.1); vec2_ <- c(0,0,1,0.3)
           mat_sex <- cbind(vec1_, vec2_)
           for (i in seq_len(n)) {
             Tp  <- exp(mat_race[, race[i] + 1] + mat_sex[, sex[i] + 1])
             Tp  <- Tp / sum(Tp)
             T[i] <- sample(T_vals, 1, prob = rep_len(Tp, k))
             Dp  <- exp(eta_); Dp <- Dp / sum(Dp)
             Delta[i] <- sample(delta_vals, 1, prob = Dp)
           }
         },
         stop("Invalid setup")
  )
  
  Indicator <- as.integer(Delta >= T)
  
  df <- data.frame(
    Delta     = Delta,
    T         = T,
    Indicator = Indicator,
    race      = race,
    sex       = sex,
    cont      = age
  )
  
  cond_pmf <- function(race, sex, cont) {
    if (!race %in% 0:2) stop("race must be 0,1,2")
    if (!sex  %in% 0:1) stop("sex must be 0,1")
    switch(as.character(setup),
           "0" = {
             pmf <- rep(1/k, k)
           },
           "1" = {
             Tp <- exp(c(2,1,0,0)); Tp <- Tp/sum(Tp)
             pmf <- rep_len(Tp, k)
           },
           "2" = {
             Dp <- exp(c(0,0.1,0.2,0.5,3,2,1,0.5,0.1))
             Dp <- Dp/sum(Dp)
             t_for_d <- pmin(4 - delta_vals + 1, 4)
             pmf <- vapply(T_vals, function(t) sum(Dp[t_for_d == t]), numeric(1))
           },
           "3" = {
             vec1  <- c(2,1,0,0); vec2 <- c(0,1,2,0); vec3 <- c(0,3,1,0)
             mat_race <- cbind(vec1, vec2, vec3)
             vec1_ <- c(0.3,1,-0.3,-0.1); vec2_ <- c(0,0,1,0.3)
             mat_sex <- cbind(vec1_, vec2_)
             Tp <- exp(mat_race[, race + 1] + mat_sex[, sex + 1])
             pmf <- rep_len(Tp / sum(Tp), k)
           }
    )
    names(pmf) <- paste0("t=", T_vals)
    pmf
  }
  
  mean_time <- function(race, sex, cont) {
    pmf <- cond_pmf(race, sex, cont)
    sum(T_vals * pmf)
  }
  
  list(
    dataframe = df,
    cond_pmf   = cond_pmf,
    mean_time  = mean_time
  )
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