library(pROC)
library(dplyr)
library(icenReg)

MIC <- function(df, method, user_formula = NULL, cov = NULL) {
  
  if(method == "np"){
    model <- ic_np(cbind(df$L, df$R))
  } else if(method %in% c("po", "ph")){
    full_formula <- as.formula(
      paste("Surv(L, R, type = 'interval2') ~", user_formula)
    )
    model <- ic_sp(full_formula, data = df, model = method)
  } else if(method %in% c("roc", "pred", "pred.adj")){
    if(method == "roc"){
      library(pROC)
      rocobj <- roc(Indicator ~ Delta, data = df, quiet = TRUE)
      cuty <- coords(rocobj, x = "best", input = "threshold", ret = "threshold", 
                     best.method = "youden", transpose = TRUE)
      threshold <- sample(cuty, 1)
      model <- list(method = "roc", threshold = threshold)
      
    } else if(method == "pred"){
      prev <- mean(df$Indicator)
      logods <- log(prev/(1 - prev))
      fit <- glm(Indicator ~ Delta, data = df, family = "binomial")
      c_ <- coef(fit)[1]
      b_ <- coef(fit)[2]
      threshold <- (logods - c_) / b_
      model <- list(method = "pred", fit = fit, threshold = threshold)
      
    } else if(method == "pred.adj"){
      prev <- mean(df$Indicator)
      logods <- log(prev/(1 - prev))
      fit <- glm(Indicator ~ Delta, data = df, family = "binomial")
      c_ <- coef(fit)[1]
      b_ <- coef(fit)[2]
      mic.pred <- (logods - c_) / b_
      SD <- sd(df$Delta)
      Cor <- cor(df$Delta, df$Indicator)
      Scf <- 0.09 * SD + 0.103 * SD * Cor
      threshold <- mic.pred - Scf * logods
      model <- list(method = "pred.adj", fit = fit, threshold = threshold)
    }
  } else {
    stop("Unsupported method: ", method)
  }
  
  # For the new methods, we use the estimated threshold as meanT.
  if(method %in% c("roc", "pred", "pred.adj")){
    meanT <- model$threshold
  } else {
    if(is.null(cov)){
      n <- nrow(df)
      v <- numeric(n)
      cov_vec <- extract_covs(user_formula)
      for (i in 1:n){
        v[i] <- mean_from_surv(model, df[i, cov_vec])
      }
      meanT <- mean(v)
    } else {
      meanT <- mean_from_surv(model, cov)
    }
  }
  
  meanD <- mean(df$Delta)
  
  result <- list(
    model = model,
    meanT = meanT,
    meanD = meanD
  )
  
  return(result)
}

tv_distance <- function(p, q) {
  0.5 * sum(abs(p - q))
}


# Assume the scores are in {1, 2, ..., score_max}

#For NPMLE 
#max_score set to 5 always



#Assume the scores are in {1, 2, ..., score_max}
transform_df <- function(df, score_max){
  
  #check if there is non-sensicle data 
  if(any(df$Delta < 1 & df$Indicator == 1)){ 
    warning("Exists at least one patient with no score improvement but reports feeling significantly better")
  }
  if(any(df$Delta == (score_max-1) & df$Indicator == 0)){ 
    warning("Exists at least one patient with max score improvement but does not feel significantly better")
  }
  
  df$L <- ifelse(df$Indicator == 1, 0, df$Delta) # (0, score_max - 1] not (1, score_max - 1]
  df$R <- ifelse(df$Indicator == 1, df$Delta, score_max-1) 
  
  return(df)
}

pmf <- function(model, cov = NULL, max_time) {
  
  model <- getSCurves(model, cov)  
  
  # Extract Turnbull intervals
  Tbull_ints <- model$Tbull_ints
  
  # Extract survival curve estimates
  S_curve <- model$S_curves[[1]]
  k <- length(S_curve)
  
  # Compute probability masses from the survival drops
  p <- numeric(k)
  p[1] <- 1 - S_curve[1]
  if (k > 1) {
    p[2:k] <- S_curve[1:(k - 1)] - S_curve[2:k]
  }
  
  # Initialize a fixed-length PMF
  pmf_fixed <- numeric(max_time)
  
  # Observed event times (excluding first element)
  observed_times <- Tbull_ints[, 2][-1]
  
  # Assign observed probabilities to correct indices
  for (i in seq_along(observed_times)) {
    time_point <- observed_times[i]
    if (time_point <= max_time && time_point >= 1) {
      pmf_fixed[time_point] <- p[i + 1] # +1 due to removal of first element earlier
    }
  }
  
  result <- list("x" = 1:max_time, "y" = pmf_fixed)
  return(result)
}


mean_from_surv <- function(model, cov = NULL) {
  
  model <- getSCurves(model, cov)  
  
  # Extract Turnbull intervals (a k x 2 matrix)
  Tbull_ints <- model$Tbull_ints
  
  # Extract survival curve estimates
  S_curve <- model$S_curves[[1]] #rather strange, but this is how it goes
  k <- length(S_curve)
  
  # Compute probability masses from the survival drops:
  p <- numeric(k)
  p[1] <- 1 - S_curve[1]
  if (k > 1) {
    p[2:k] <- S_curve[1:(k - 1)] - S_curve[2:k]
  }
  
  points <- Tbull_ints[, 2]
  mean_event_time <- sum(p * points)
  
  return(mean_event_time)
  
}

#written by ChatGPT
#Utility function that extracts the cov from user_formula
#Used for computing E[T] =E[E[T|X]]
extract_covs <- function(formula_str) {
  # Remove spaces
  formula_str <- gsub(" ", "", formula_str)
  
  # Split the string by '+'
  terms <- unlist(strsplit(formula_str, "\\+"))
  
  # For each term, remove the 'factor()' wrapper if it exists
  terms <- sapply(terms, function(x) {
    if (grepl("^factor\\(", x)) {
      # Extract the content inside factor()
      sub("^factor\\((.*)\\)$", "\\1", x)
    } else {
      x
    }
  })
  return(terms)
}
