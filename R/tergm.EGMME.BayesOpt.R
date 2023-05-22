# Function outline:

# Take in a ergm state and other stuff as in EGMME.GD
# Get an Estiamte of the network params (this should be done in GD so can copy)
# Calculate the mahalanobis distance/cost function (this should be done in GD so can copy)
# Repeat

tergm.EGMME.bayesOpt <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  
  theta <- theta0
  
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accomodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  # always collect if monitoring model is passed
  control$collect <- TRUE
  control$changes <- FALSE
  
  control$time.burnin <- 1
  control$time.interval <- 1
  control$time.samplesize <- 1

  
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(model$etamap$etalength), model.mon$nw.stats-model.mon$target.stats))
  
  # call many times
  # re-run every time update theta
  eta <- ergm.eta(theta, model$etamap)
  eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
  # browser()
  z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
  
  #sometging like this... z.stats$matrix
  # stats should already be centred so no need to do subtraction
  
  # cost function takes theta, map onto eta, call mcmc slave, get statistics, get mahanob distance
}