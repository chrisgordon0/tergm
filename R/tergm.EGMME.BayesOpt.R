# Function outline:

# Take in a ergm state and other stuff as in EGMME.GD
# Get an Estiamte of the network params (this should be done in GD so can copy)
# Calculate the mahalanobis distance/cost function (this should be done in GD so can copy)
# Repeat

tergm.EGMME.GD <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){

states <- replicate(nthreads(control),
            {
              list(nw=nw,
                   eta = ergm.eta(theta0, model$etamap),
                   nw.diff  = nw.diff)
            },
            simplify=FALSE
)


out <- parallel::clusterApply(ergm.getCluster(control), seq_along(states), function(i) tergm_MCMC_sample(states[[i]]$nw, model, model.mon, proposal, eta=states[[i]]$eta, control=control.phase1, verbose=verbose))

print(out)


}