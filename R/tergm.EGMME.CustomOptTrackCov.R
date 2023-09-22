tergm.EGMME.customOptTrackCov <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){

  starting_thetas <- numeric()
  final_thetas <- numeric()
  final_distances <- numeric()
  
  plot(COTV.create_even_grid(10, list(c(-10,0), c(-3,5))))
  COTV.create_even_grid(10, list(c(-10,0), c(-3,5)))
  for (k in 1:30) {
  
  set.seed(k)
  theta0 <- unname(theta0)
  theta0 <- c(rnorm(1,-2.94439, sd=3), rnorm(1, 1, sd=3))
  starting_thetas <- append(starting_thetas, theta0)
  theta_bounds <- list(c(theta0[1]-1, theta0[1]+1), c(theta0[2]-1, theta0[2]+1))
  print(theta_bounds)
  

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
  #theta_bounds <- list(c(-4, -2), c(0, 2))
  
  # theta0 <- unname(theta0)
  # theta_bounds <- list(c(theta0[1]-1, theta0[1]+1), c(theta0[2]-1, theta0[2]+1))

  threshold <- 0.5
  
  # Initial Run
  des <- COTV.getDesign(8, theta_bounds)
  thetas <- as.matrix(des, ncol=2)
  target_stats <- COTV.getTargetStatSample(des, model, model.mon, current_ergm_state, control, verbose)
  thetas <- thetas[rep(1:nrow(thetas), each = 2), ]
  target_stats_GPs <- COTV.fitGauProcToTargetStats(target_stats, thetas)
  cov_matrix_star <- COTV.getCovMatrixFromTargetStats(theta0, model, model.mon, current_ergm_state, control, verbose)
  obj_fn <- COTV.predictObjectiveWrapper(target_stats_GPs, theta_bounds, 2, cov_matrix_star, search_min=FALSE)
  
  print("Objective Function:")
  print(obj_fn)
  
  min_dist <- min(obj_fn[,3])
  min_theta <- obj_fn[which.min(obj_fn[,3]),1:2]
  
  print("Min dist")
  print(min_dist)
  print("Min theta")
  print(min_theta)
  
  min_theta_history <- obj_fn[which.min(obj_fn[,3]),]
  
  i <- 0
  while(i < 6 && min_dist > 0.1) {
    
    result <- COTV.checkIfThetaCloseToBoundsAndReturnNewBounds(theta_bounds, min_theta, threshold)
    
    theta_bounds <- result$new_bounds_total
    new_area <- result$new_area
    added_bounds_1 <- result$added_bounds_1
    added_bounds_2 <- result$added_bounds_2
    

    new_thetas_to_sample <- NULL
    if (added_bounds_1) {
      area <- list(new_area[[1]], theta_bounds[[2]])
      print("area1")
      print(area)
      # search these bounds in the design matrix thing
      des <- COTV.getDesign(8, area)
      new_thetas_to_sample <- as.matrix(des, ncol=2)
    }
    if (added_bounds_2) {
      area <- list(theta_bounds[[1]], new_area[[2]])
      print("area2")
      print(area)
      # search these bounds in the design matrix thing
      des <- COTV.getDesign(8, area)
      if(!is.null(new_thetas_to_sample)){
        new_thetas_to_sample <- rbind(new_thetas_to_sample, as.matrix(des, ncol=2))
      }
      else {
        new_thetas_to_sample <- as.matrix(des, ncol=2)
      }
    }
    
    print("New theta bounds")
    print(theta_bounds)
    
    new_theta <- COTV.findNextThetaToSample(obj_fn)
    new_theta <- matrix(new_theta, nrow=1, ncol=2)
    new_thetas_to_sample <- rbind(new_thetas_to_sample, new_theta)
    
    print("New thetas to sample")
    print(new_thetas_to_sample)
    
    target_stats_new <- COTV.getTargetStatSample(new_thetas_to_sample, model, model.mon, current_ergm_state, control, verbose)
    
    target_stats <- rbind(target_stats, target_stats_new)
    
    new_thetas_to_sample <- new_thetas_to_sample[rep(1:nrow(new_thetas_to_sample), each = 2), ]
    thetas <- rbind(thetas, new_thetas_to_sample)
    
    target_stats_GPs <- COTV.fitGauProcToTargetStats(target_stats, thetas)
    obj_fn <- COTV.predictObjectiveWrapper(target_stats_GPs, theta_bounds, 2, cov_matrix_star, search_min=FALSE)
    
    print("New Obj Fn")
    print(obj_fn)
    
    min_dist <- min(obj_fn[,3])
    min_theta <- obj_fn[which.min(obj_fn[,3]),1:2]
    
    print("Min dist")
    print(min_dist)
    print("Min theta")
    print(min_theta)
    
    min_theta_history <- rbind(min_theta_history,  obj_fn[which.min(obj_fn[,3]),])
    
    print(min_theta_history)
    
    if (!(added_bounds_1 || added_bounds_2)) {
      prev_min_theta <- min_theta_history[nrow(min_theta_history),1:2]
      theta_dist <- sqrt((min_theta[1]-prev_min_theta[1])^2+(min_theta[2]-prev_min_theta[2])^2)
      print(theta_dist)
      if (theta_dist < 0.2) {
        print("Min dist was too close, searching the next minimum")
        combinations <- expand.grid(theta_bounds)
        result <- obj_fn[obj_fn[,1] %in% combinations$Var1 & obj_fn[,2] %in% combinations$Var2, 3]
        
        min_index <- which.min(result)
        min_theta <- c(combinations$Var1[min_index], combinations$Var2[min_index])
        print("new min theta")
        print(min_theta)
      }
    }
    
    i <- i+1
  }
  print("END")
  print(min_theta)
  print(min_dist)
  final_thetas <- append(final_thetas, min_theta)
  final_distances <- append(final_distances, min_dist)
  } # end of big for loop
  
  print("ALL RESULTS")
  print("Starting Thetas")
  print(starting_thetas)
  print("Final Thetas")
  print(final_thetas)
  print("Final Distances")
  print(final_distances)
  
  res <- cbind(matrix(starting_thetas, ncol=2, byrow=TRUE), matrix(final_thetas, ncol=2, byrow=TRUE), final_distances)
  print(res)
  write.csv(res, file='/Users/chris/Documents/uni_work/tergm_data/customOptResults/results.csv', row.names=FALSE)
  return()

}

COTV.checkIfThetaCloseToBoundsAndReturnNewBounds <- function(theta_bounds, min_theta, threshold) {
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



COTV.fitGauProcToTargetStats <- function(target_stats, sampled_thetas) {
  models <- list()
  for (i in 1:ncol(target_stats)) {
    #gp <- list(gpr(target_stats[,i], sampled_thetas, gamma=1.5))
    gp <- list(gpr(target_stats[,i], sampled_thetas, Cov="matern", nu=3/2, mean=1))
    models <- append(models, gp)
  }
  return(models)
}


COTV.predictObjectiveWrapper <- function(target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star, search_min=FALSE) {
  result <- COTV.predictObjective(40, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star)
  
  if (!search_min) {
    return(result)
  }
  
  obj_fn <- result
  
  min_theta <- obj_fn[which.min(obj_fn[,3]),1:2]
  theta_bounds <- list(c(min_theta[1]-0.25, min_theta[1]+0.25), c(min_theta[2]-0.25, min_theta[2]+0.25))
    
  result <- COTV.predictObjective(10, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star)
  
  result <- rbind(obj_fn, result)
  
  return(result)
}

COTV.predictObjective <- function(points_per_dim, target_stats_GPs, theta_bounds, n_target_stats, cov_matrix_star) {
  theta_grid <- COTV.create_even_grid(points_per_dim, theta_bounds)
  
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
    dist <- COTV.mahalanobisDist(cov_matrix, target_stats_predictions[i, 1:n_target_stats])
    distances <- c(distances, dist)
    
    var <- 10^10
    if (dist != 10^10) {
      var <- COTV.getVarianceOfDistance(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats, cov_matrix_star)
    }
    variances <- c(variances, var)
    
    #COTV.printIfNegative(cov_matrix, target_stats_predictions, i, n_target_stats)
  }
  
  result <- cbind(theta_grid, distances, variances)
  
  #write.csv(result, file='temp.csv', row.names=FALSE)
  
  return(result)
}

COTV.getVarianceOfDistance <- function(cov_matrix, target_stats_predictions, target_stats_pred_sd, i, n_target_stats, cov_matrix_star) {
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
  
  #corr_mat <- COTV.getCorrelationMatrix(target_stats_predictions[i,], target_stats_pred_sd[i,], n_target_stats)

  diag(cov_matrix_star) <- target_stats_pred_sd[i,]^2
  #print(cov_matrix_star)
  #cov_matrix_star <- diag(target_stats_pred_sd[i,]^2)
  #print(cov_matrix_star)
  var_dist <- t(grad_dist) %*% cov_matrix_star %*% grad_dist
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
  
  # print("ATTACH EI")
  # print(cbind(objective_fn, all_ei))
  # print("CHOSEN")
  # print(max_improvement)
  # print(best_theta)
  
  return(best_theta)
}


# For minimization
COTV.EI <- function(mu, sigma, ybest) {
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
COTV.getTargetStatSample <- function(theta_values_to_sample, model, model.mon, current_ergm_state, control, verbose) {
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
    
    var1 <- (statsmatrix[,1] - mean(statsmatrix[,1]))^2
    mat <- cbind(mat, matrix(c(mean(var1) - sqrt(var(var1)/2), mean(var1) + sqrt(var(var1)/2)), ncol=1, nrow=2))
    
    var2 <- (statsmatrix[,2] - mean(statsmatrix[,2]))^2
    mat <- cbind(mat, matrix(c(mean(var2) - sqrt(var(var2)/2), mean(var2) + sqrt(var(var2)/2)), ncol=1, nrow=2))
    
    cov1 <- (statsmatrix[,1]-mean(statsmatrix[,1]))*(statsmatrix[,2]-mean(statsmatrix[,2]))
    mat <- cbind(mat, matrix(c(mean(cov1)-sqrt(var(cov1)/2), mean(cov1)+sqrt(var(cov1)/2)), ncol=1, nrow=2))
    mat <- cbind(mat, matrix(c(mean(cov1)-sqrt(var(cov1)/2), mean(cov1)+sqrt(var(cov1)/2)), ncol=1, nrow=2))
    
    # print(statsmatrix)
    # print(var(statsmatrix[,1]))
    # print(mat)
    # print("---------------------")
    return(mat)
  }
  
  stopCluster(cl)
  registerDoSEQ()

  return(target_stats)
}

COTV.getCovMatrixFromTargetStats <- function(theta, model, model.mon, current_ergm_state, control, verbose) {

  statsmatrix <- COTV.getMCMCSample(theta, model, model.mon, current_ergm_state, control, verbose)
  
  var1 <- (statsmatrix[,1]-mean(statsmatrix[,1]))^2
  var2 <- (statsmatrix[,2]-mean(statsmatrix[,2]))^2
  cov1 <- (statsmatrix[,1]-mean(statsmatrix[,1]))*(statsmatrix[,2]-mean(statsmatrix[,2]))
  
  return(cov(cbind(statsmatrix, var1, var2, cov1, cov1)))
}


COTV.getMCMCSample <- function(theta, model, model.mon, current_ergm_state, control, verbose) {
  set.seed(1)
  
  eta <- ergm::ergm.eta(theta, model$etamap)
  eta.comb <- c(eta, numeric(model.mon$etamap$etalength))
  
  z <- tergm::tergm_MCMC_slave(current_ergm_state, eta.comb, control, verbose)
  statsmatrix <- z$statsmatrix[,-seq_len(nparam(model, canonical = TRUE))]
  
  return(statsmatrix)
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
