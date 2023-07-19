tergm.EGMME.customOpt <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){

  #print(-seq_len(nparam(model, canonical = TRUE)))
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  control <- set_control_options(control)
  
  theta_bounds <- list(c(-10, 10), c(-10, 10))
  #theta_bounds <- list(c(-6, -5), c(2, 3))
  #theta_bounds <- list(c(-10,10))
  points_per_dim <- 5
  
  theta_values_to_sample <- create_even_grid(points_per_dim, theta_bounds)
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                  stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  num_target_stats <- 2 # (I think this is nparam(model, canonical = TRUE))
  
  # For fitting it with the GP that is fit to the variance and covariance
  merged_target_stats <- generateInitialTargetStatSample_trackCovariance(theta_values_to_sample, points_per_dim, model, model.mon, current_ergm_state, control, verbose)
  
  library(plot3D)
  
  x <- theta_values_to_sample[,1]
  y <- theta_values_to_sample[,2]
  z <- merged_target_stats[,2]
  
  surf3D(x = x, y = y, z = z, colkey=TRUE, bty="b2",
         phi = 40, theta = 30)
  return()
  
  print(cbind(theta_values_to_sample, merged_target_stats))
  target_stats_GPs <- fitGauProcToTargetStats(merged_target_stats, theta_values_to_sample)
  predictObjective_trackCovariance(target_stats_GPs, theta_bounds, num_target_stats)
  return()
  
  # For fitting it the normal way
  #merged_target_stats <- generateInitialTargetStatSample(theta_values_to_sample, points_per_dim, model, model.mon, current_ergm_state, control, verbose)
  
  # Normal way but load from file
  #merged_target_stats <- getMergedTargetStatsFromFile(theta_bounds)
  #print(merged_target_stats)
  #theta_values_to_sample <- merged_target_stats[,1]
  #merged_target_stats <- merged_target_stats[,-1]
  #print(merged_target_stats)
  
  # Not needed
  #merged_target_stats <- cbind(theta_values_to_sample, merged_target_stats)
  
  #write.table(merged_target_stats, file="./R/target_stats_new.csv", row.names = F, col.names = F)
  
  
  #print("THETAS")
  #print(theta_values_to_sample)

  #print("TARGET STATS, CROSS TARGET STATS")
  #print(merged_target_stats)
  
  #target_stats_GPs <- fitGauProcToTargetStatsParallel(merged_target_stats, theta_values_to_sample)
  #target_stats_GPs <- fitGauProcToTargetStats(merged_target_stats, theta_values_to_sample)
  target_stats_GPs <- fitMultiGauProcToTargetStats(merged_target_stats, theta_values_to_sample)
  return()
  
  # Plot the target stats with their GP approximation
  # This is all just plot output for debugging
  n_plots <- ncol(merged_target_stats)
  layout(matrix(1:n_plots, nrow = 2, byrow = TRUE))

  for (i in 1:ncol(merged_target_stats)) {
    gp <- target_stats_GPs[[i]]
    start = theta_bounds[[1]][1]
    end = theta_bounds[[1]][2]
    grid <- seq(start, end, by = 0.01)
    ytest <- gprPredict(train=gp, inputNew=grid)$pred.mean
    par(mar = c(4, 4, 2, 1))
    plot(grid, ytest, type = "l", ylim = range(ytest, merged_target_stats[, i]),
         xlab = "Grid", ylab = "Values", col = "red")
    lines(theta_values_to_sample, merged_target_stats[, i], type = 'l', col = 'blue')
    #lines(theta_values_to_sample, merged_target_stats[, i], type = 'p', col = 'blue')

    title(main = paste("Plot", i))
  }

  #layout(1)
  
  #objective_fn <- predictObjective(target_stats_GPs, theta_bounds, num_target_stats)
  
  #print(objective_fn)
  #next_theta <- findNextThetaToSample(objective_fn, theta_bounds)
  

  
  return(NULL)
}

generateInitialTargetStatSample_trackCovariance <- function(theta_values_to_sample, points_per_dim, model, model.mon, current_ergm_state, control, verbose) {
  num_cores <- 4
  cl <- makeCluster(num_cores, outfile="")
  registerDoParallel(cl)
  
  target_stats <- vector("list", nrow(theta_values_to_sample))
  
  target_stats <- foreach(i = 1:nrow(theta_values_to_sample), .combine = "rbind") %dopar% {
    set.seed(1)
    theta <- theta_values_to_sample[i,]
    eta <- ergm.eta(theta, model$etamap)
    eta.comb <- c(eta, numeric(model.mon$etamap$etalength))

    z <- tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)

    statsmatrix <- z$statsmatrix[,-seq_len(nparam(model, canonical = TRUE))]

    variances <- c(mean(statsmatrix[,1]^2) - mean(statsmatrix[,1])^2, mean(statsmatrix[,2]^2) - mean(statsmatrix[,2])^2)
    covariances <- c(mean(statsmatrix[,1]*statsmatrix[,2])-mean(statsmatrix[,1])*mean(statsmatrix[,2]), mean(statsmatrix[,2]*statsmatrix[,1])-mean(statsmatrix[,2])*mean(statsmatrix[,1]))

    output <- c(colMeans(statsmatrix), variances, covariances)
    output <- matrix(output, nrow=1)
    return(output)
  }
  stopCluster(cl)
  registerDoSEQ()
  
  return(target_stats)
}


# Note: cross_target_stats has columns in order
# 1x1, 1x2, 1x3, 2x1, 2x2, 2x3, etc
# where 1,2,3 are the columns of target_stats
# Merged Target Stats is like
# 1, 2, 3, 1x1, 1x2, 1x3, 2x1, 2x2, 2x3, etc
generateInitialTargetStatSample <- function(theta_values_to_sample, points_per_dim, model, model.mon, current_ergm_state, control, verbose) {
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
    statsmatrix <- z$statsmatrix[,-seq_len(nparam(model, canonical = TRUE))]
    output <- c(colMeans(statsmatrix), getCrossTargetStatistics(statsmatrix))
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
    gp <- list(gpr(target_stats[,i], sampled_thetas, gamma=1.5)) #kernel="" I think uses RBF kernel by default
    models <- append(models, gp)
  }
  return(models)
}

fitGauProcToTargetStatsParallel <- function(target_stats, sampled_thetas) {
  num_cores <- 4
  cl <- makeCluster(num_cores, outfile="")
  clusterExport(cl=cl, c("gpr"))
  registerDoParallel(cl)
  
  models <- foreach (i = 1:ncol(target_stats), .combine='c', .multicombine=TRUE) %dopar% {
    gp <- list(gpr(target_stats[,i], sampled_thetas, gamma=2)) #kernel="" I think uses RBF kernel by default
    return(gp)
  }
  
  return(models)
}

fitMultiGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  target_stats <- matrix(unlist(target_stats), ncol=length(target_stats))

  replicated_sampled_thetas <- replicate(6, sampled_thetas)
  
  input_list <- lapply(1:ncol(replicated_sampled_thetas), function(i) replicated_sampled_thetas[,i])
  
  data <- list(input=input_list, response=lapply(1:ncol(target_stats), function(i) matrix(target_stats[,i], ncol=1)))

  # print(lapply(1:nrow(target_stats), function(i) target_stats[i,]))
  
  m <- mgpr(Data = data)
  
  input1 <- -5
  data_new <- list()
  data_new$input <- list(input1)
  mgprPredict(m, DataNew=data_new)

  #print(target_stats)
}

 predictObjective_trackCovariance <- function(target_stats_GPs, theta_bounds, num_target_stats) {
  theta_grid <- create_even_grid(20, theta_bounds)

  target_stats_predictions <- NULL
  for (i in 1:length(target_stats_GPs)) {
    target_stats_predictions <- cbind(target_stats_predictions, gprPredict(train=target_stats_GPs[[i]], inputNew=theta_grid)$pred.mean)
  }

  distances <- numeric()

  min_distance <- Inf
  min_distance_theta <- NULL
  idx_of_min <- 0
  for (i in 1:nrow(target_stats_predictions)) {
    cov_matrix <- matrix(c(target_stats_predictions[i, 3], target_stats_predictions[i,5], target_stats_predictions[i,6], target_stats_predictions[i,4]), ncol=2, nrow=2, byrow=TRUE)
    #print(cov_matrix)
    
    distances <- c(distances, mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]))
    
    if (mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]) < min_distance) {
      min_distance <- mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats])
      min_distance_theta <- theta_grid[i,]
      idx_of_min <- i
    }
    
  }
  #print(distances)
  print("MIN DISTANCE")
  print(min_distance)
  print("MIN DISTANCE THETA")
  print(min_distance_theta)
  print("COV MATRIX")
  print(matrix(c(target_stats_predictions[idx_of_min, 3], target_stats_predictions[idx_of_min,5], target_stats_predictions[idx_of_min,6], target_stats_predictions[idx_of_min,4]), ncol=2, nrow=2, byrow=TRUE))
  # 
  # plot(theta_grid[1:150], distances[1:150])
  # lines(theta_grid, distances, type='l', col='blue')
  # 
  # plot(theta_grid, distances)
  # lines(theta_grid, distances, type='l', col='blue')
  
  return(cbind(theta_grid, distances))
}

predictObjective <- function(target_stats_GPs, theta_bounds, num_target_stats) {
  theta_grid <- create_even_grid(500, theta_bounds)
  #print(theta_grid)
  
  target_stats_predictions <- NULL
  for (i in 1:length(target_stats_GPs)) {
    target_stats_predictions <- cbind(target_stats_predictions, gprPredict(train=target_stats_GPs[[i]], inputNew=theta_grid)$pred.mean)
  }
  # print(target_stats_predictions)
  
  distances <- numeric()
  
  min_distance <- Inf
  min_distance_theta <- NULL
  
  for (i in 1:nrow(target_stats_predictions)) {
    cov_matrix <- getCovarianceMatrix(target_stats_predictions[i,], num_target_stats)
    
    if (theta_grid[i,] >= -9.6 && theta_grid[i,] <= -9.4) {
      print(theta_grid[i,])
      print(mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]))
      print(target_stats_predictions[i,])
      print(cov_matrix)
    }
    distances <- c(distances, mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]))
    #print(distances)
    #return()
    
    if (mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats]) < min_distance) {
      min_distance <- mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:num_target_stats])
      min_distance_theta <- theta_grid[i]
      
      # print("New Min:")
      # print(min_distance)
      # print(min_distance_theta)
    }
    
  }
  # print("ALL DISTANCES")
  chosen_dist <- Inf
  chosen_theta_row <- -1
  # for (i in 1:nrow(theta_grid)) {
  #   if (chosen_dist > distances[i] && distances[i] >= 0) {
  #     chosen_dist <- distances[i]
  #     chosen_theta_row <- i
  #   }
  #   print(theta_grid[i,])
  #   print(distances[i])
  # }
  # print("CHOSEN MIN")
  # print(chosen_dist)
  # print(theta_grid[chosen_theta_row,])
  # print(min(distances))
  # print(theta_grid[match(min(distances), distances),1])
  # print(theta_grid[match(min(distances), distances),2])
  # print(theta_grid)
  # print(distances)
  plot(theta_grid[1:300], distances[1:300])
  lines(theta_grid, distances, type='l', col='blue')

  plot(theta_grid, distances)
  lines(theta_grid, distances, type='l', col='blue')
  return(cbind(theta_grid, distances))
}

# Requires expected value of: first moment, second moment and cross moments
# target_stats: a vector (1 row) containing the above
mahalanobisDist <- function(cov_matrix, target_statistics) {
  # print(cov_matrix)
  
  # check for singular matrix
  if (rcond(cov_matrix) < .Machine$double.eps) {
    return(10^3)
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

# work in progress
# plot the GP at some point to ensure its right
findNextThetaToSample <- function(objective_fn, theta_bounds) {
  # Squish the objective
  # Need to do this as otherwise the GP for the objective isn't accurate (needs to inverse a matrix and makes floating point errors)
  #shift <- abs(min(objective_fn[,ncol(objective_fn)])) + 1
  #objective_fn[,ncol(objective_fn)] <- log(objective_fn[,ncol(objective_fn)] + shift)
  #objective_fn[,ncol(objective_fn)] <- matrix(scale(objective_fn[,ncol(objective_fn)]), ncol=1)
  
  objective_gp <- gpr(objective_fn[,ncol(objective_fn)], objective_fn[,-ncol(objective_fn)]) #kernel="" I think uses RBF kernel by default
  
  ybest <- min(objective_fn[,ncol(objective_fn)])
  possible_thetas <- create_even_grid(50, theta_bounds)
  
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
  plot(possible_thetas, gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.mean, type='l', col='black', ylab="GP Prediction for Objective (standardised)")
  lines(possible_thetas, gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.mean+1.96*gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.sd, col='blue')
  lines(possible_thetas, gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.mean-1.96*gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.sd, col='blue')
  #print(gprPredict(train=objective_gp, inputNew=possible_thetas)$pred.sd)
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




######## HELPERS ########

set_control_options <- function(control) {
  control$collect <- TRUE # always collect if monitoring model is passed
  control$changes <- FALSE
  control$time.burnin <- 100
  control$time.interval <- 10
  control$time.samplesize <- 200
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

getMergedTargetStatsFromFile <- function(theta_bounds) {
   merged_target_stats <- read.csv("./R/target_stats_new.csv", header=FALSE, sep=" ")
   subset_rows <- merged_target_stats[merged_target_stats[,1] >= theta_bounds[[1]][1] & merged_target_stats[,1] <= theta_bounds[[1]][2], ]
   return(subset_rows)
}