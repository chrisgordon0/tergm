tergm.EGMME.customOptTrackCovNew <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){
  
  
  theta0 <- unname(theta0)
  theta_bounds <- list(c(theta0[1]-1, theta0[1]+1),
                       c(theta0[2]-1, theta0[2]+1),
                       c(theta0[3]-1, theta0[3]+1),
                       c(theta0[4]-1, theta0[4]+1),
                       c(theta0[5]-1, theta0[5]+1))
  n_target_stats <- 5
  
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  current_ergm_state <- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(model$etamap$etalength), model.mon$nw.stats - model.mon$target.stats))
  
  control <- COTVN.set_control_options(control)
  
  threshold <- 0.5
  
  # Initial Run
  des <- COTVN.getDesign(20, theta_bounds)
  thetas <- as.matrix(des, ncol=2)
  
  
  target_stats <- COTVN.getTargetStatSample(des, model, model.mon, current_ergm_state, control, verbose)
  thetas <- thetas[rep(1:nrow(thetas), each = 2), ]
  target_stats_GPs <- COTVN.fitGauProcToTargetStats(target_stats, thetas)
  
  cov_matrix_star <- COTVN.getCovMatrixFromTargetStats(theta0, model, model.mon, current_ergm_state, control, verbose)
  obj_fn <- COTVN.predictObjectiveWrapper(target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star, search_min=FALSE)
  
  print("Objective Function:")
  print(obj_fn)
  
  min_dist <- min(obj_fn[,3])
  min_theta <- obj_fn[which.min(obj_fn[,3]),1:5]
  
  print("Min dist")
  print(min_dist)
  print("Min theta")
  print(min_theta)
  
  return()
  
  min_theta_history <- obj_fn[which.min(obj_fn[,3]),]
  
  i <- 0
  while(i < 6 && min_dist > 0.1) {
    
    result <- COTVN.checkIfThetaCloseToBoundsAndReturnNewBounds(theta_bounds, min_theta, threshold)
    
    theta_bounds <- result$new_bounds_total
    new_area <- result$new_area
    added_bounds_1 <- result$added_bounds_1
    added_bounds_2 <- result$added_bounds_2
    
    
    new_thetas_to_sample <- NULL
    if (added_bounds_1) {
      area <- list(new_area[[1]], theta_bounds[[2]])
      ##print("area1")
      #print(area)
      # search these bounds in the design matrix thing
      des <- COTVN.getDesign(8, area)
      new_thetas_to_sample <- as.matrix(des, ncol=2)
    }
    if (added_bounds_2) {
      area <- list(theta_bounds[[1]], new_area[[2]])
      #print("area2")
      #print(area)
      # search these bounds in the design matrix thing
      des <- COTVN.getDesign(8, area)
      if(!is.null(new_thetas_to_sample)){
        new_thetas_to_sample <- rbind(new_thetas_to_sample, as.matrix(des, ncol=2))
      }
      else {
        new_thetas_to_sample <- as.matrix(des, ncol=2)
      }
    }
    
    
    new_theta <- COTVN.findNextThetaToSample(obj_fn)
    new_theta <- matrix(new_theta, nrow=1, ncol=2)
    new_thetas_to_sample <- rbind(new_thetas_to_sample, new_theta, min_theta)
    
    
    target_stats_new <- COTVN.getTargetStatSample(new_thetas_to_sample, model, model.mon, current_ergm_state, control, verbose)
    
    target_stats <- rbind(target_stats, target_stats_new)
    
    new_thetas_to_sample <- new_thetas_to_sample[rep(1:nrow(new_thetas_to_sample), each = 2), ]
    thetas <- rbind(thetas, new_thetas_to_sample)
    
    target_stats_GPs <- COTVN.fitGauProcToTargetStats(target_stats, thetas)
    obj_fn <- COTVN.predictObjectiveWrapper(target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star, search_min=FALSE)
    
    
    min_dist <- min(obj_fn[,3])
    min_theta <- obj_fn[which.min(obj_fn[,3]),1:2]
    
    
    min_theta_history <- rbind(min_theta_history,  obj_fn[which.min(obj_fn[,3]),])

    if (!(added_bounds_1 || added_bounds_2)) {
      prev_min_theta <- min_theta_history[nrow(min_theta_history),1:2]
      theta_dist <- sqrt((min_theta[1]-prev_min_theta[1])^2+(min_theta[2]-prev_min_theta[2])^2)

        if (theta_dist < 0.2) {
        combinations <- expand.grid(theta_bounds)
        result <- obj_fn[obj_fn[,1] %in% combinations$Var1 & obj_fn[,2] %in% combinations$Var2, 3]
        
        min_index <- which.min(result)
        min_theta <- c(combinations$Var1[min_index], combinations$Var2[min_index])

      }
    }
    
    i <- i+1
  }
  
  return(min_theta)
  
}

COTVN.checkIfThetaCloseToBoundsAndReturnNewBounds <- function(theta_bounds, min_theta, threshold) {
  added_bounds_1 <- FALSE
  added_bounds_2 <- FALSE
  alt_theta_bounds <- theta_bounds
  
  
  if (min_theta[1] < (theta_bounds[[1]][1] + threshold)) {
    alt_theta_bounds[[1]][2] <- theta_bounds[[1]][1]
    theta_bounds[[1]][1] <- theta_bounds[[1]][1] - 1
    alt_theta_bounds[[1]][1] <- theta_bounds[[1]][1]
    added_bounds_1 <- TRUE
  }
  
  if (min_theta[1] > (theta_bounds[[1]][2] - threshold)) {
    alt_theta_bounds[[1]][1] <- theta_bounds[[1]][2]
    theta_bounds[[1]][2] <- theta_bounds[[1]][2] + 1
    alt_theta_bounds[[1]][2] <- theta_bounds[[1]][2]
    added_bounds_1 <- TRUE
  }
  
  if (min_theta[2] < (theta_bounds[[2]][1] + threshold)) {
    alt_theta_bounds[[2]][2] <- theta_bounds[[2]][1]
    theta_bounds[[2]][1] <- theta_bounds[[2]][1] - 1
    alt_theta_bounds[[2]][1] <- theta_bounds[[2]][1]
    added_bounds_2 <- TRUE
  }
  
  if (min_theta[2] > (theta_bounds[[2]][2] - threshold)) {
    alt_theta_bounds[[2]][1] <- theta_bounds[[2]][2]
    theta_bounds[[2]][2] <- theta_bounds[[2]][2] + 1
    alt_theta_bounds[[2]][2] <- theta_bounds[[2]][2]
    added_bounds_2 <- TRUE
  }
  
  return(list(new_bounds_total = theta_bounds, new_area = alt_theta_bounds, added_bounds_1 = added_bounds_1, added_bounds_2 = added_bounds_2))
}

# For now its ok but later I can only fit a subset and copy the dupes
COTVN.fitGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  models <- list()
  for (i in 1:ncol(target_stats)) {
    gp <- list(gpr(target_stats[,i], sampled_thetas, Cov="matern", nu=3/2, mean=1))
    models <- append(models, gp)
  }
  return(models)
}


COTVN.predictObjectiveWrapper <- function(target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star, search_min=FALSE) {
  result <- COTVN.predictObjective(5, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star)
  
  if (!search_min) {
    return(result)
  }
  
  obj_fn <- result
  
  min_theta <- obj_fn[which.min(obj_fn[,3]),1:2]
  theta_bounds <- list(c(min_theta[1]-0.25, min_theta[1]+0.25),
                       c(min_theta[2]-0.25, min_theta[2]+0.25),
                       c(min_theta[3]-0.25, min_theta[3]+0.25),
                       c(min_theta[4]-0.25, min_theta[4]+0.25),
                       c(min_theta[5]-0.25, min_theta[5]+0.25))
  
  result <- COTVN.predictObjective(10, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star)
  
  result <- rbind(obj_fn, result)
  
  return(result)
}

COTVN.predictObjective <- function(points_per_dim, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star) {
  theta_grid <- COTVN.create_even_grid(points_per_dim, theta_bounds)
  
  # Order: x_1, x_2, v_1, c_12, v_2, c_21 (where c_12 and c_21 should be equal)
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
    cov_matrix <-  matrix(target_stats_predictions[i,6:ncol(target_stats_predictions)], ncol=5, byrow=TRUE)
    dist <- COTVN.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats])
    distances <- c(distances, dist)
    
    var <- 10^10
    if (dist != 10^10) {
      var <- COTVN.getVarianceOfDistance(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats, cov_matrix_star)
    }
    variances <- c(variances, var)
    
    #COTVN.#printIfNegative(cov_matrix, target_stats_predictions, i, n_target_stats)
  }
  
  result <- cbind(theta_grid, distances, variances)
  
  #write.csv(result, file='temp.csv', row.names=FALSE)
  
  return(result)
}

COTVN.getVarianceOfDistance <- function(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats, cov_matrix_star) {
  X <- target_stats_predictions[i, 1:5]
  sigma <- cov_matrix
  sigma_inv <- solve(sigma)
  
  X <- as.matrix(X, ncol=1, nrow=5)
  
  grad_X <- 2 * t(X) %*% sigma_inv
  grad_Sigma <- -sigma_inv %*% X %*% t(X) %*% sigma_inv
  
  
  grad_dist <- matrix(grad_X, nrow=5)
  grad_dist <- rbind(grad_dist, matrix(grad_Sigma, ncol=1))
  
  
  #corr_mat <- COTVN.getCorrelationMatrix(target_stats_predictions[i,], target_stats_pred_sd[i,], n_target_stats)
  
  diag(cov_matrix_star) <- target_stats_pred_sd[i,]^2
  ##print(cov_matrix_star)
  #cov_matrix_star <- diag(target_stats_pred_sd[i,]^2)
  ##print(cov_matrix_star)
  var_dist <- t(grad_dist) %*% cov_matrix_star %*% grad_dist
  return(var_dist)
}

COTVN.getCorrelationMatrix <- function(target_stats_predictions, target_stats_pred_sd, n_target_stats) {
  corr_mat <- diag(target_stats_pred_sd^2)
  
  for (i in 1:length(target_stats_predictions)) {
    for (j in 1:length(target_stats_predictions)) {
      if (i == j) {
        next
      }
      corr_mat[i,j] <- COTVN.getCov(target_stats_predictions[i], target_stats_predictions[j], target_stats_pred_sd[i], target_stats_pred_sd[j]) 
      corr_mat[j,i] <- corr_mat[i,j]
    }
  }
  
  return(corr_mat)
}

COTVN.getCov <- function(x_1, x_2, sd_1, sd_2) {
  x_1_sample <- rnorm(10000, mean=x_1, sd=sd_1)
  x_2_sample <- rnorm(10000, mean=x_2, sd=sd_2)
  return(cov(x_1_sample,x_2_sample))
}

COTVN.printIfNegative <- function(cov_matrix, target_stats_predictions, i, n_target_stats) {
  if (COTVN.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats]) < 0) {
    #print(COTVN.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats]))
    #print(cov_matrix)
    #print("")
  }
}

COTVN.findNextThetaToSample <- function(objective_fn) {
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
    ei <- COTVN.EI(pred, sqrt(var), ybest)
    if (ei > max_improvement) {
      max_improvement <- ei
      best_theta <- current_theta
    }
    all_ei <- c(all_ei, ei)
  }
  
  # #print("ATTACH EI")
  # #print(cbind(objective_fn, all_ei))
  # #print("CHOSEN")
  # #print(max_improvement)
  # #print(best_theta)
  
  return(best_theta)
}


# For minimization
COTVN.EI <- function(mu, sigma, ybest) {
  # for now
  if (is.na(sigma)) {
    return(10^10)
  }
  if (sigma == 0) {
    return(0)
  }
  gamma <- (ybest - mu) / sigma
  phi <- pnorm(gamma)
  return(sigma * (gamma * phi + dnorm(gamma)))
}


# Calculates the target stats and their variance/covariances. 
# Also finds the two representative points.
COTVN.getTargetStatSample <- function(theta_values_to_sample, model, model.mon, current_ergm_state, control, verbose) {
  num_cores <- 2
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
                            
                            # This should make it v1, c12, c13, ..., c1n, v2, c21, c22, ..., c2n, ...
                            for(i in 1:ncol(statsmatrix)) {
                              for (j in 1:ncol(statsmatrix)) {
                                if (i == j) {
                                  variance <- (statsmatrix[,i] - mean(statsmatrix[,i]))^2
                                  mat <- cbind(mat, matrix(c(mean(variance) - sqrt(var(variance)/2), mean(variance) + sqrt(var(variance)/2)), ncol=1, nrow=2))
                                }
                                else {
                                  covariance <- (statsmatrix[,i]-mean(statsmatrix[,i]))*(statsmatrix[,j]-mean(statsmatrix[,j]))
                                  mat <- cbind(mat, matrix(c(mean(covariance) - sqrt(var(covariance)/2), mean(covariance) + sqrt(var(covariance)/2)), ncol=1, nrow=2))
                                }
                              }
                            }
                            
                            return(mat)
                          }
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(target_stats)
}

COTVN.getCovMatrixFromTargetStats <- function(theta, model, model.mon, current_ergm_state, control, verbose) {
  
  statsmatrix <- COTVN.getMCMCSample(theta, model, model.mon, current_ergm_state, control, verbose)
  final_mat <- statsmatrix
  for(i in 1:ncol(statsmatrix)) {
    for (j in 1:ncol(statsmatrix)) {
      if (i == j) {
        variance <- (statsmatrix[,i] - mean(statsmatrix[,i]))^2
        final_mat <- cbind(final_mat, variance)
      }
      else {
        covariance <- (statsmatrix[,i]-mean(statsmatrix[,i]))*(statsmatrix[,j]-mean(statsmatrix[,j]))
        final_mat <- cbind(final_mat, covariance)
      }
    }
  }
  return(cov(final_mat))
}


COTVN.getMCMCSample <- function(theta, model, model.mon, current_ergm_state, control, verbose) {
  set.seed(1)
  
  eta <- ergm::ergm.eta(theta, model$etamap)
  eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
  
  z <- tergm::tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
  statsmatrix <- z$statsmatrix[,-seq_len(nparam(model, canonical = TRUE))]
  
  return(statsmatrix)
}


COTVN.mahalanobisDist <- function(cov_matrix, target_statistics) {
  # check for singular matrix
  if (rcond(cov_matrix) <= .Machine$double.eps) {
    return(10^10)
  }
  
  inv_cov_matrix <- solve(cov_matrix)
  dist <- mahalanobis(target_statistics, rep(0, length(target_statistics)), inv_cov_matrix, inverted=TRUE)
  
  if (dist < 0) {
    # MAYBE TRY NEARPD as SNEARPD sometimes the square root is undefined.
    # cov_matrix <- statnet.common::snearPD(cov_matrix)$mat
    tryCatch(
      {
        cov_matrix <- Matrix::nearPD(cov_matrix)$mat
        cov_matrix <- as.matrix(cov_matrix)
        
        inv_cov_matrix <- solve(cov_matrix)
        dist <- mahalanobis(target_statistics, rep(0, length(target_statistics)), inv_cov_matrix, inverted=TRUE)
      },
      error = function(e) {
        dist <- 10^10
      }
    )
  }
  
  if (dist < 0) {return (10^10)}
  
  return(dist)
}


# Helper to clean up main function
COTVN.getDesign <- function(n, theta_bounds) {
  params = makeParamSet(
    makeNumericParam("theta1", lower=theta_bounds[[1]][1], upper=theta_bounds[[1]][2]),
    makeNumericParam("theta2", lower=theta_bounds[[2]][1], upper=theta_bounds[[2]][2]),
    makeNumericParam("theta3", lower=theta_bounds[[3]][1], upper=theta_bounds[[3]][2]),
    makeNumericParam("theta4", lower=theta_bounds[[4]][1], upper=theta_bounds[[4]][2]),
    makeNumericParam("theta5", lower=theta_bounds[[5]][1], upper=theta_bounds[[5]][2])
  )
  
  return(ParamHelpers::generateDesign(n = n, par.set = params, fun = lhs::randomLHS))
}


COTVN.create_even_grid <- function(points_per_dim, bounds) {
  seqs <- list()
  for (i in 1:length(bounds)) {
    seqs <- append(seqs, list(seq(bounds[[i]][1], bounds[[i]][2], length.out = points_per_dim)))
  }
  grid <- expand.grid(seqs)
  return(matrix(unlist(grid), ncol=ncol(grid)))
}


COTVN.set_control_options <- function(control) {
  control$collect <- TRUE # always collect if monitoring model is passed
  control$changes <- FALSE
  control$time.burnin <- 1000
  control$time.interval <- 1
  control$time.samplesize <- 100
  return(control)
}
