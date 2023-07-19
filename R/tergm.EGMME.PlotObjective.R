tergm.EGMME.plotObjective <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  
  
  #distances <- read.csv("./R/objective_distances_new.csv", header=FALSE, sep=" ")
  #plot(distances[,1], distances[,2], xlim=c(-10,0), ylim=c(0,1000))

  #all_objective_distances <- read.csv("./R/objective_distances_predicted_1unit.csv", header=FALSE, sep=" ")
  #lines(all_objective_distances[,1], all_objective_distances[,2], col='blue', type='p')
  #lines(all_objective_distances[,1], all_objective_distances[,2], col='blue')


  #return()
  
  
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  control <- set_control_options(control)
  
  theta_values <- c(-9.473684) #seq(-10, 10, by=0.025)
  #theta_values <- c(-5.979899)#, 2.5)
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  
  set.seed(1)
  theta <- theta_values[1]
  eta <- ergm.eta(theta, model$etamap)
  eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
  z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
  
  target_statistics <- z$statsmatrix[,-(1:2)]
  print(target_statistics)
  cov_matrix <- cov(target_statistics)
  # check for singular matrix
  if (rcond(cov_matrix) < .Machine$double.eps) {
    return(10^6)
  }
  print(cov_matrix)
  inv_cov_matrix <- solve(cov_matrix)
  
  dist <- mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix)
  print(dist)  
  return()
  # distances <- numeric()
  
  # for (i in 1:length(theta_values)) {
  #   distances <- c(distances, optimCostFunction1(theta_values[i], control, verbose, model, model.mon))
  # }
  # 
  # plot(theta_values, distances)
  # lines(theta_values, distances, type='l')
  # 
  # distances <- matrix(distances, ncol=1)
  # distances <- cbind(theta_values, distances)
  
  # write.table(distances, file="./R/objective_distances.csv", row.names = F, col.names = F)
  
  # i=1
  # theta <- theta_values[i]
  # eta <- ergm.eta(theta, model$etamap)
  # eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
  # z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
  # #print(z$statsmatrix[,-(1:2)])
  # #print(z$statsmatrix)
  # 
  # target_statistics <- z$statsmatrix[,-(1:2)]
  # cov_matrix <- cov(target_statistics)
  # #print(cov_matrix)
  # # check for singular matrix
  # if (rcond(cov_matrix) < .Machine$double.eps) {
  #   return(10^13)
  # }
  # 
  # inv_cov_matrix <- solve(cov_matrix)
  # 
  # dist <- mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix)
  # 
  # output <- matrix(dist, nrow=1)
  # print(output)
  # return()

  
  num_cores <- 4
  cl <- makeCluster(num_cores, outfile="")
  clusterExport(cl=cl, c("getCrossTargetStatistics"))
  registerDoParallel(cl)
  
  distances <- vector("list", length(theta_values))
  
  distances <- foreach(i = 1:length(theta_values), .combine = "rbind") %dopar% {
    set.seed(1)
    theta <- theta_values[i]
    eta <- ergm.eta(theta, model$etamap)
    eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
    z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
    
    target_statistics <- z$statsmatrix[,-(1:2)]
    cov_matrix <- cov(target_statistics)
    # check for singular matrix
    if (rcond(cov_matrix) < .Machine$double.eps) {
      return(10^6)
    }
    
    inv_cov_matrix <- solve(cov_matrix)
    
    dist <- mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix)

    output <- matrix(dist, nrow=1)
    return(output)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  distances <- cbind(theta_values, distances)
  write.table(distances, file="./R/objective_distances_new.csv", row.names = F, col.names = F)
  
  plot(distances[,1], distances[,2])
  
  return(NULL)
}