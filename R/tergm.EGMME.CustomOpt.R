
tergm.EGMME.customOpt <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  

  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  control <- set_control_options(control)
  
  #theta_bounds <- list(c(-6, -4), c(1,3))
  theta_bounds <- list(c(-6, -3))
  points_per_dim <- 20
  theta_values_to_sample <- create_even_grid(points_per_dim, theta_bounds)

  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                  stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  # Figure out how to get this earlier
  num_target_stats <- 3
  
  merged_target_stats <- generateInitialTargetStatSample(theta_values_to_sample, points_per_dim, theta0, model, model.mon, current_ergm_state, control, verbose)

  # # Plot the target stats against their index
  # # Set the layout for the grid
  # n_plots <- ncol(merged_target_stats)
  # layout(matrix(1:n_plots, nrow = 2, byrow = TRUE))
  # 
  # # Loop through each column of merged_target_stats
  # for (i in 1:ncol(merged_target_stats)) {
  #   # Set the plot area for the current plot
  #   plot_index <- i
  #   par(mar = c(4, 4, 2, 1))
  #   plot(seq(1, nrow(merged_target_stats), by = 1), merged_target_stats[, i],
  #        xlab = "Index", ylab = paste("Column", i))
  #   
  #   # Add a title for each plot
  #   title(main = paste("Plot", i))
  # }
  
  # Reset the layout
  layout(1)
  
  
  print("THETAS")
  print(theta_values_to_sample)

  print("TARGET STATS, CROSS TARGET STATS")
  print(merged_target_stats)

  target_stats_GPs <- fitGauProcToTargetStats(merged_target_stats, theta_values_to_sample)

  objective_fn <- predictObjective(target_stats_GPs, theta_bounds, num_target_stats)
  next_theta <- findNextThetaToSample(objective_fn, theta_bounds)
  
  # Plot the target stats with their GP approximation
  # This is all just plot output for debugging
  n_plots <- ncol(merged_target_stats)
  layout(matrix(1:n_plots, nrow = 2, byrow = TRUE))

  for (i in 1:ncol(merged_target_stats)) {
    gp <- target_stats_GPs[[i]]
    grid <- seq(-6, -3, by = 0.01)
    ytest <- gprPredict(train=gp, inputNew=grid)$pred.mean
    par(mar = c(4, 4, 2, 1))
    plot(grid, ytest, type = "l", ylim = range(ytest, merged_target_stats[, i]),
         xlab = "Grid", ylab = "Values", col = "red")
    lines(theta_values_to_sample, merged_target_stats[, i], type = 'l', col = 'blue')
    lines(theta_values_to_sample, merged_target_stats[, i], type = 'p', col = 'blue')

    title(main = paste("Plot", i))
  }

  layout(1)
  
  return(NULL)
}


# Note: cross_target_stats has columns in order
# 1x1, 1x2, 1x3, 2x1, 2x2, 2x3, etc
# where 1,2,3 are the columns of target_stats
# Merged Target Stats is like
# 1, 2, 3, 1x1, 1x2, 1x3, 2x1, 2x2, 2x3, etc
generateInitialTargetStatSample <- function(theta_values_to_sample, points_per_dim, theta0, model, model.mon, current_ergm_state, control, verbose) {
  num_cores <- 4
  cl <- makeCluster(num_cores, outfile="")
  clusterExport(cl=cl, c("getCrossTargetStatistics"))
  registerDoParallel(cl)

  target_stats <- vector("list", nrow(theta_values_to_sample))
  
  target_stats <- foreach(i = 1:nrow(theta_values_to_sample), .combine = "rbind") %dopar% {
    set.seed(1)
    theta <- theta_values_to_sample[i,]
    eta <- ergm.eta(theta, model$etamap)
    eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
    z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
    output <- c(colMeans(z$statsmatrix), getCrossTargetStatistics(z$statsmatrix))
    output <- matrix(output, nrow=1)
    return(output)
  }
  print(target_stats)
  stopCluster(cl)
  registerDoSEQ()

  return(target_stats)
}

fitGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  models <- list()
  for (i in 1:ncol(target_stats)) {
    gp <- list(gpr(target_stats[,i], sampled_thetas)) #kernel="" I think uses RBF kernel by default
    models <- append(models, gp)
  }
  return(models)
}

# work in progress
# plot the GP at some point to ensure its right
findNextThetaToSample <- function(objective_fn, theta_bounds) {
  # Squish the objective
  # Need to do this as otherwise the GP for the objective isn't accurate (needs to inverse a matrix and makes floating point errors)
  #shift <- abs(min(objective_fn[,ncol(objective_fn)])) + 1
  #objective_fn[,ncol(objective_fn)] <- log(objective_fn[,ncol(objective_fn)] + shift)
  objective_fn[,ncol(objective_fn)] <- matrix(scale(objective_fn[,ncol(objective_fn)]), ncol=1)
  
  objective_gp <- gpr(objective_fn[,ncol(objective_fn)], objective_fn[,-ncol(objective_fn)]) #kernel="" I think uses RBF kernel by default
  ybest <- min(objective_fn[,ncol(objective_fn)])
  possible_thetas <- create_even_grid(1000, theta_bounds)

  best_theta <- NULL
  max_improvement <- -Inf
  all_ei <- numeric()
  for (i in 1:nrow(possible_thetas)) {
    current_theta <- possible_thetas[i,]
    prediction <- gprPredict(train=objective_gp, inputNew=current_theta)
    ei <- expectedImprovement(prediction$pred.mean, prediction$pred.sd, ybest)
    if (ei > max_improvement) {
      max_improvement <- ei
      best_theta <- current_theta
    }
    all_ei <- c(all_ei, ei)
  }
  gprPredict(train=objective_gp, inputNew=possible_thetas)
  # Plot the GP
  plot(possible_thetas, gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.mean, type='l', col='black', ylab="GP Prediction for Objective (normalized objective)")
  # Plot the EI
  plot(possible_thetas, all_ei, type='l', col='blue', ylab="Expected Improvement")
  print("CHOSEN")
  print(max_improvement)
  print(best_theta)
  return(best_theta)
}

# For minimization
expectedImprovement <- function(mu, sigma, ybest) {
  #print(mu)
  #print(sigma)
  #print(ybest)
  # Assuming NaN is 0 but need to check this because a point against itself
  # should still have uncertainty
  #if (is.na(sigma)) {
  #  return(0)
  #}
  if (sigma == 0) {
    return(0)
  }
  gamma <- (ybest - mu) / sigma
  phi <- pnorm(gamma)
  return(sigma * (gamma * phi + dnorm(gamma)))
}


predictObjective <- function(target_stats_GPs, theta_bounds, num_target_stats) {
  theta_grid <- create_even_grid(20, theta_bounds)
  
  target_stats_predictions <- NULL
  for (i in 1:length(target_stats_GPs)) {
    target_stats_predictions <- cbind(target_stats_predictions, gprPredict(train=target_stats_GPs[[i]], inputNew=theta_grid)$pred.mean)
  }
  print(target_stats_predictions)
  
  distances <- numeric()  
  
  min_distance <- Inf
  min_distance_theta <- NULL
  
  for (i in 1:nrow(target_stats_predictions)) {
    cov_matrix <- getCovarianceMatrix(target_stats_predictions[i,], num_target_stats)
    distances <- c(distances, mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]))
    
    if (mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]) < min_distance) {
      min_distance <- mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats])
      min_distance_theta <- theta_grid[i]
      
      # print("New Min:")
      # print(min_distance)
      # print(min_distance_theta)
    }
    
  }

  # print(distances)
  plot(theta_grid, distances)
  lines(theta_grid, distances, type='l', col='blue')
  
  return(cbind(theta_grid, distances))
}

# Requires expected value of: first moment, second moment and cross moments
# target_stats: a vector (1 row) containing the above
mahalanobisDist <- function(cov_matrix, target_statistics) {
  
  # check for singular matrix
  if (rcond(cov_matrix) < .Machine$double.eps) {
    return(10^13)
  }

  inv_cov_matrix <- solve(cov_matrix)
  dist <- mahalanobis(target_statistics, rep(0, length(target_statistics)), inv_cov_matrix)
  return(dist)
}

# Requires expected value of: first moment, second moment and cross moments
# target_stats: a vector (1 row) containing the above
getCovarianceMatrix <- function(target_stats, num_target_stats) {
  cov_matrix <- matrix(0, nrow = num_target_stats, ncol = num_target_stats)
  
  for (i in 1:num_target_stats) {
    for (j in 1:num_target_stats) {
      if (i == j) {
        # E(t(Y)^2) - E(t(Y))^2
        cov_matrix[i,i] <- target_stats[i*(num_target_stats+1)] - target_stats[i]^2
      } else {
        # E(t_i(Y)t_j(Y))- E(t_i(Y))E(t_j(Y))
        cov_matrix[i,j] <- target_stats[i*num_target_stats + j] - target_stats[i]*target_stats[j]
      }
    }
  }
  return(cov_matrix)
}

getCrossTargetStatistics <- function(target_statistics, num_target_stats) {
  cross_target_stats <- numeric()
  n_cols <- ncol(target_statistics)
  column_combinations <- expand.grid(1:n_cols, 1:n_cols)
  for (i in 1:nrow(column_combinations)) {
    col1 <- column_combinations[i, 2]
    col2 <- column_combinations[i, 1]
    
    cross_target_stats <- c(cross_target_stats, mean(target_statistics[,col1] * target_statistics[,col2]))
  }
  return(cross_target_stats)
}

######## HELPERS ########

set_control_options <- function(control) {
  control$collect <- TRUE # always collect if monitoring model is passed
  control$changes <- FALSE
  control$time.burnin <- 10
  control$time.interval <- 10
  control$time.samplesize <- 100
  return(control)
}

create_even_grid <- function(points_per_dim, bounds) {
  seqs <- list()
  for (i in 1:length(bounds)) {
    seqs <- append(seqs, list(seq(bounds[[i]][1], bounds[[i]][2], length.out = points_per_dim)))
  }
  grid <- expand.grid(seqs)
  return(matrix(unlist(grid), ncol=ncol(grid)))
}