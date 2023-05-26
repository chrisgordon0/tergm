# Function outline:

# Take in a ergm state and other stuff as in EGMME.GD
# Get an Estiamte of the network params (this should be done in GD so can copy)
# Calculate the mahalanobis distance/cost function (this should be done in GD so can copy)
# Repeat

tergm.EGMME.bayesOpt <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){

  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal

  # always collect if monitoring model is passed
  control$collect <- TRUE
  control$changes <- FALSE
  
  control$time.burnin <- 1
  control$time.interval <- 1
  control$time.samplesize <- 100

  
  # Must assign globals here, and after the optimization I should clean them up
  # and remove them from the environment.
  # The reason I must make them globals is because the cost function for
  # Bayesian Optimization can only have one parameter.
  global.control <<- control
  global.theta <<- theta0
  global.nw <<- nw
  global.model <<- model
  global.model.comb <<- model.comb
  global.model.mon <<- model.mon
  global.proposal <<- proposal
  global.verbose <<- verbose
  
  # Need to do this as the cost function cannot be a vector
  # (limitation in the rBayesianOptimization package)
  # TODO: Maybe change which optimization package we use.
  bounds <- create_bounds(theta0)

  optim_result <- BayesianOptimization(optimCostFunctionWrapper,
                                       bounds = bounds,
                                       init_points = 2, n_iter = 2,
                                       acq = "ucb", kappa = 1, eps = 0.0,
                                       verbose = TRUE)
  
  optimal_inputs <- optim_result$Best_Par
  optimal_value <- optim_result$Best_Value
  
  print(optimal_inputs)
  print(optimal_value)
    
  # At the end do something like this to clean up:
  # rm(global.control) ... etc
  
  
}

create_bounds <- function(theta) {
  bounds <- list()
  bounds_length <- length(theta)
  
  for (i in 1:bounds_length) {
    bounds[[paste0("x", i)]] <- c(-1, 1)
  }
  
  return(bounds)
}

mahalanobisDist <- function(target_statistics) {
  
  cov_matrix <- cov(target_statistics)

  if (rcond(cov_matrix) < 0.01) {
    return(10000)
  }
  

  cov <<- cov_matrix
  
  inv_cov_matrix <- solve(cov_matrix)
  
  inv <<- inv_cov_matrix
  
  dist <- sqrt(mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix))
  return(dist)
}

optimCostFunction <- function(theta) {
  
  current_ergm_state <- ergm_state(global.nw, model=global.model.comb, proposal=global.proposal,
                                   stats=c(numeric(global.model$etamap$etalength), global.model.mon$nw.stats - global.model.mon$target.stats))
  
  # call many times
  # re-run every time update theta
  eta <- ergm.eta(theta, global.model$etamap)
  eta.comb <- c(eta, numeric(global.model.mon$etamap$etalength))

  z <- tergm_MCMC_slave(current_ergm_state, eta.comb, global.control, global.verbose)
  
  # optional to change ergm state every time
  current_ergm_state <- z$state
  
  target_statistics <- z$statsmatrix
  
  targets <<- target_statistics
  
  # For debugging
  mahalanob <<- mahalanobisDist(target_statistics)
  
  return(list(Score = mahalanobisDist(target_statistics), Pred = 0.1))  
}

optimCostFunctionWrapper <- function(...) {
  theta <- c(...)
  return(optimCostFunction(theta))
}

