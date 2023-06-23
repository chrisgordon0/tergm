
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
  
  target_stats <- numeric()
  # Note: cross_target_stats has columns in order
  # 1x1, 1x2, 1x3, 2x1, 2x2, 2x3, etc
  # where 1,2,3 are the columns of target_stats
  cross_target_stats <- numeric()
  
  # can remove this in the future
  sampled_thetas <- numeric()
  
  # How do I get this before hand?
  num_target_stats <- 3
  
  for (i in 1:nrow(theta_values_to_sample)) {
    print("done")
    theta <- theta_values_to_sample[i,]
    sampled_thetas <- c(sampled_thetas, theta)

    eta <- ergm.eta(theta, model$etamap)
    eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
    z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)

    target_stats <- c(target_stats, colMeans(z$statsmatrix))
    cross_target_stats <- c(cross_target_stats, getCrossTargetStatistics(z$statsmatrix))
  }
  
  sampled_thetas <- matrix(sampled_thetas, ncol=length(theta0), byrow=TRUE) # check byrow
  target_stats <- matrix(target_stats, nrow=points_per_dim^length(theta0), byrow=TRUE)
  cross_target_stats <- matrix(cross_target_stats, nrow=points_per_dim^length(theta0), byrow=TRUE)
  
  merged_target_stats <- cbind(target_stats, cross_target_stats)

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
  print(sampled_thetas)

  print("TARGET STATS, CROSS TARGET STATS")
  print(merged_target_stats)

  target_stats_GPs <- fitGauProcToTargetStats(merged_target_stats, sampled_thetas)

  predictObjective(target_stats_GPs, theta_bounds, num_target_stats)

  # This is all just plot output for debugging
  n_plots <- ncol(merged_target_stats)
  layout(matrix(1:n_plots, nrow = 2, byrow = TRUE))

  for (i in 1:ncol(merged_target_stats)) {
    gp <- target_stats_GPs[[i]]
    grid <- seq(-6, -4, by = 0.01)
    ytest <- kernlab::predict(gp, grid)
    par(mar = c(4, 4, 2, 1))
    plot(grid, ytest, type = "l", ylim = range(ytest, merged_target_stats[, i]),
         xlab = "Grid", ylab = "Values", col = "red")
    lines(sampled_thetas, merged_target_stats[, i], type = 'l', col = 'blue')
    lines(sampled_thetas, merged_target_stats[, i], type = 'p', col = 'blue')

    title(main = paste("Plot", i))
  }

  layout(1)
  
  
  return(NULL)
}


fitGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  models <- list()
  for (i in 1:ncol(target_stats)) {
    gp <- gausspr(sampled_thetas, target_stats[,i], var=1, kernel=) # I think uses RBF kernel by default
    models <- append(models, gp)
  }
  return(models)
}


predictObjective <- function(target_stats_GPs, theta_bounds, num_target_stats) {
  theta_grid <- create_even_grid(100, theta_bounds)
  
  print(theta_grid)
  
  target_stats_predictions <- NULL
  for (i in 1:length(target_stats_GPs)) {
    target_stats_predictions <- cbind(target_stats_predictions, kernlab::predict(target_stats_GPs[[i]], theta_grid))
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
      
      print("New Min:")
      print(min_distance)
      print(min_distance_theta)
    }
    
  }

  print(distances)
  plot(theta_grid, distances)
  lines(theta_grid, distances, type='l', col='blue')
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