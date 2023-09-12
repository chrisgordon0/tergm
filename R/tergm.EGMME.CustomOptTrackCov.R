tergm.EGMME.customOptTrackCov <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  control <- COTV.set_control_options(control)
  
  # Note theta0 is (-2.9444,1.000) for SimpleTest
  
  #theta_bounds <- list(c(-6, -2), c(0, 4))
  #theta_bounds <- list(c(-10,0), c(0,10))
  theta_bounds <- list(c(-5, -1), c(-1, 3))
  

  des <- COTV.getDesign(15, theta_bounds)
  target_stats <- COTV.getTargetStatSample(des, model, model.mon, current_ergm_state, control, verbose)

  
  thetas <- des[rep(1:nrow(des), each = 2), ]
  
  print(cbind(thetas,target_stats))

  target_stats_GPs <- COTV.fitGauProcToTargetStats(target_stats, thetas)

  obj_fn <- COTV.predictObjective(target_stats_GPs, theta_bounds, 2)

  COTV.findNextThetaToSample(obj_fn)
  
  return()

}

COTV.fitGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  models <- list()
  for (i in 1:ncol(target_stats)) {
    #gp <- list(gpr(target_stats[,i], sampled_thetas, gamma=1.5))
    gp <- list(gpr(target_stats[,i], sampled_thetas, Cov="matern", nu=3/2))
    models <- append(models, gp)
  }
  return(models)
}


COTV.predictObjective <- function(target_stats_GPs, theta_bounds, n_target_stats) {
  theta_grid <- create_even_grid(10, theta_bounds)
  
  # Order: x_1, x_2, v_1, v_2, c_12, c_21 (where c_12 and c_21 should be equal)
  target_stats_predictions <- NULL
  target_stats_pred_sd <- NULL
  for (i in 1:length(target_stats_GPs)) {
    prediction <- gprPredict(train=target_stats_GPs[[i]], inputNew=theta_grid)
    target_stats_predictions <- cbind(target_stats_predictions, prediction$pred.mean)
    target_stats_pred_sd <- cbind(target_stats_pred_sd, prediction$pred.sd)
    
  }

  
  distances <- numeric()
  variances <- numeric()
  for (i in 1:nrow(target_stats_predictions)) {
    cov_matrix <- matrix(c(target_stats_predictions[i, 3], target_stats_predictions[i,5], target_stats_predictions[i,6], target_stats_predictions[i,4]), ncol=2, nrow=2, byrow=TRUE)
    print(cov_matrix)
    dist <- COTV.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats])
    distances <- c(distances, dist)
    
    var <- 10^10
    if (dist != 10^10) {
      var <- COTV.getVarianceOfDistance(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats)
    }
    variances <- c(variances, var)
    
    #COTV.printIfNegative(cov_matrix, target_stats_predictions, i, n_target_stats)
  }
  
  result <- cbind(theta_grid, distances, variances)
  
  #write.csv(result, file='temp.csv', row.names=FALSE)
  
  return(result)
}

COTV.getVarianceOfDistance <- function(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats) {
  X <- target_stats_predictions[i, 1:2]
  sigma <- cov_matrix
  sigma_inv <- solve(sigma)
  
  X <- as.matrix(X, ncol=1, nrow=2)

  grad_X <- 2 * t(X) %*% sigma_inv
  grad_Sigma <- -sigma_inv %*% X %*% t(X) %*% sigma_inv
  
  
  grad_dist <- matrix(grad_X, nrow=2)
  grad_dist <- rbind(grad_dist, matrix(grad_Sigma, ncol=1))
  
  # Since the order of target_stats_predictions and target_stats_pred_sd is
  # x_1, x_2, v_1, v_2, c_12, c_21 (where c_12 and c_21 should be equal)
  # we will swap the final row and the fourth row of grad_dist as we want to make
  # grad_dist: x_1, x_2, v_1, c_12, c_21, v_2
  # into
  # grad_dist: x_1, x_2, v_1, v_2, c_12, c_21
  
  desired_order <- c(1, 2, 3, 6, 4, 5)
  grad_dist <- grad_dist[desired_order, , drop = FALSE]
  
  corr_mat <- COTV.getCorrelationMatrix(target_stats_predictions[i,], target_stats_pred_sd[i,], n_target_stats)

  var_dist <- t(grad_dist) %*% corr_mat %*% grad_dist
  return(var_dist)
}

COTV.getCorrelationMatrix <- function(target_stats_predictions, target_stats_pred_sd, n_target_stats) {
  corr_mat <- diag(target_stats_pred_sd^2)

  for (i in 1:length(target_stats_predictions)) {
    for (j in 1:length(target_stats_predictions)) {
      if (i == j) {
        next
      }
      corr_mat[i,j] <- COTV.getCov(target_stats_predictions[i], target_stats_predictions[j], target_stats_pred_sd[i], target_stats_pred_sd[j]) 
      corr_mat[j,i] <- corr_mat[i,j]
    }
  }
  
  return(corr_mat)
}

COTV.getCov <- function(x_1, x_2, sd_1, sd_2) {
  x_1_sample <- rnorm(10000, mean=x_1, sd=sd_1)
  x_2_sample <- rnorm(10000, mean=x_2, sd=sd_2)
  return(cov(x_1_sample,x_2_sample))
}

COTV.printIfNegative <- function(cov_matrix, target_stats_predictions, i, n_target_stats) {
  if (COTV.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats]) < 0) {
    print(COTV.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats]))
    print(cov_matrix)
    print("")
  }
}

COTV.findNextThetaToSample <- function(objective_fn) {
  thetas <- objective_fn[,1:2]
  distances <- objective_fn[,3]
  variances <- objective_fn[,4]

  ybest <- min(distances)

  best_theta <- NULL
  max_improvement <- -Inf
  all_ei <- numeric()
  for (i in 1:nrow(objective_fn)) {
    current_theta <- thetas[i,]
    pred <- distances[i]
    var <- variances[i]
    ei <- COTV.EI(pred, sqrt(var), ybest)
    if (ei > max_improvement) {
      max_improvement <- ei
      best_theta <- current_theta
    }
    all_ei <- c(all_ei, ei)
  }
  
  print("ATTACH EI")
  print(cbind(objective_fn, all_ei))
  print("CHOSEN")
  print(max_improvement)
  print(best_theta)
  
  return(best_theta)
}


# For minimization
COTV.EI <- function(mu, sigma, ybest) {
  if (sigma == 0) {
    return(0)
  }
  gamma <- (ybest - mu) / sigma
  phi <- pnorm(gamma)
  return(sigma * (gamma * phi + dnorm(gamma)))
}


# Calculates the target stats and their variance/covariances. 
# Also finds the two representative points.
COTV.getTargetStatSample <- function(theta_values_to_sample, model, model.mon, current_ergm_state, control, verbose) {
  num_cores <- 4
  cl <- makeCluster(num_cores, outfile="")
  registerDoParallel(cl)
  
  target_stats <- vector("list", nrow(theta_values_to_sample)*2)
  
  target_stats <- foreach(i = 1:nrow(theta_values_to_sample), .combine = "rbind",
                        .packages=c("ergm", "statnet.common", "tergm")) %dopar% {
                          
    set.seed(1)
                          
    theta <- theta_values_to_sample[i,]
    eta <- ergm::ergm.eta(theta, model$etamap)
    eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
    
    z <- tergm::tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
    statsmatrix <- z$statsmatrix[,-seq_len(nparam(model, canonical = TRUE))]
    mat <- NULL
    for (i in 1:ncol(statsmatrix)) {
      sample <- matrix(c(mean(statsmatrix[,i]) - sqrt(var(statsmatrix[,i])/2),
                            mean(statsmatrix[,i]) + sqrt(var(statsmatrix[,i])/2)),
                          ncol=1, nrow=2)
      if(is.null(mat)) {
        mat <- sample
      }
      else {
        mat <- cbind(mat, sample)
      }
    }
    
    var1 <- (statsmatrix[,1] - mean(statsmatrix[,1]))^2
    mat <- cbind(mat, matrix(c(mean(var1) - sqrt(var(var1)/2), mean(var1) + sqrt(var(var1)/2)), ncol=1, nrow=2))
    
    var2 <- (statsmatrix[,2] - mean(statsmatrix[,2]))^2
    mat <- cbind(mat, matrix(c(mean(var2) - sqrt(var(var2)/2), mean(var2) + sqrt(var(var2)/2)), ncol=1, nrow=2))
    
    cov1 <- (statsmatrix[,1]-mean(statsmatrix[,1]))*(statsmatrix[,2]-mean(statsmatrix[,2]))
    mat <- cbind(mat, matrix(c(mean(cov1)-sqrt(var(cov1)/2), mean(cov1)+sqrt(var(cov1)/2)), ncol=1, nrow=2))
    mat <- cbind(mat, matrix(c(mean(cov1)-sqrt(var(cov1)/2), mean(cov1)+sqrt(var(cov1)/2)), ncol=1, nrow=2))
    
    return(mat)
  }
  
  stopCluster(cl)
  registerDoSEQ()

  return(target_stats)
}


COTV.mahalanobisDist <- function(cov_matrix, target_statistics) {
  # check for singular matrix
  if (rcond(cov_matrix) <= .Machine$double.eps) {
    return(10^10)
  }

  inv_cov_matrix <- solve(cov_matrix)
  dist <- mahalanobis(target_statistics, rep(0, length(target_statistics)), inv_cov_matrix, inverted=TRUE)
  
  if (dist < 0) {
    # MAYBE TRY NEARPD as SNEARPD sometimes the square root is undefined.
    #cov_matrix <- statnet.common::snearPD(cov_matrix)$mat
    cov_matrix <- Matrix::nearPD(cov_matrix)$mat
    cov_matrix <- as.matrix(cov_matrix)
    
    inv_cov_matrix <- solve(cov_matrix)
    dist <- mahalanobis(target_statistics, rep(0, length(target_statistics)), inv_cov_matrix, inverted=TRUE)
  }

  return(dist)
}


# Helper to clean up main function
COTV.getDesign <- function(n, theta_bounds) {
  params = makeParamSet(
    makeNumericParam("theta1", lower=theta_bounds[[1]][1], upper=theta_bounds[[1]][2]),
    makeNumericParam("theta2", lower=theta_bounds[[2]][1], upper=theta_bounds[[2]][2])
  )
  
  return(ParamHelpers::generateDesign(n = n, par.set = params, fun = lhs::randomLHS))
}


COTV.create_even_grid <- function(points_per_dim, bounds) {
  seqs <- list()
  for (i in 1:length(bounds)) {
    seqs <- append(seqs, list(seq(bounds[[i]][1], bounds[[i]][2], length.out = points_per_dim)))
  }
  grid <- expand.grid(seqs)
  return(matrix(unlist(grid), ncol=ncol(grid)))
}


COTV.set_control_options <- function(control) {
  control$collect <- TRUE # always collect if monitoring model is passed
  control$changes <- FALSE
  control$time.burnin <- 200
  control$time.interval <- 10
  control$time.samplesize <- 200
  return(control)
}
