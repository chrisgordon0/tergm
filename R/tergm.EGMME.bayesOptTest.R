# Function outline:

# Take in a ergm state and other stuff as in EGMME.GD
# Get an Estiamte of the network params (this should be done in GD so can copy)
# Calculate the mahalanobis distance/cost function (this should be done in GD so can copy)
# Repeat

tergm.EGMME.bayesOptTest <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  # always collect if monitoring model is passed
  control$collect <- TRUE
  control$changes <- FALSE
  
  control$time.burnin <- 1000
  control$time.interval <- 1000
  control$time.samplesize <- 10000
  
  #plotObjectiveWithTheta(theta0, nw, model, model.mon, control, proposal, verbose)
  #return()
  
  #theta <- c(-5.188136, 2.114250)
  #theta <- c(-4.810679, 1.890441)
  #theta <- c(-5.320798, 2.551421)
  theta <- c(-5.207, 2.208)
  dist <- optimCostFunction(theta)
  
  #print("Inside Test")
  print(dist)
  #print("Outside Test")
  return(dist)
}



plotObjectiveWithTheta <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE) {
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  
  theta1 <- seq(from=-6, to=6, by=0.5)
  theta2 <- seq(from=-6, to=6, by=0.5)
  res = matrix(nrow=length(theta1), ncol=length(theta2))
  print(res)
  i_ = 1
  for (i in theta1) {
    j_ = 1
    for (j in theta2) {
      theta <- c(i,j)
      # call many times
      # re-run every time update theta
      eta <- ergm.eta(theta, model$etamap)
      eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
      
      z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
      
      # optional to change ergm state every time (doesn't change anything)
      # current_ergm_state <- z$state
      
      target_statistics <- z$statsmatrix
      
      
      # For debugging
      mahalanob <- mahalanobisDist(target_statistics)
      
      cov_matrix <- cov(target_statistics)
      
      # check for singular matrix
      if (rcond(cov_matrix) < .Machine$double.eps) {
        print("Computationally Signular")
        print(10^13)
        res[i,j] = 10^13
        next
      }
      
      inv_cov_matrix <- solve(cov_matrix)
      
      dist <- sqrt(mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix))
      res[i,j] = dist
      print(dist)
    }
  }
  res <<- res
  print(res)
}