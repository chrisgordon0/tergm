# Function outline:

# Take in a ergm state and other stuff as in EGMME.GD
# Get an Estiamte of the network params (this should be done in GD so can copy)
# Calculate the mahalanobis distance/cost function (this should be done in GD so can copy)
# Repeat

tergm.EGMME.bayesOpt <- function(theta0, nw, model, model.mon, control, proposal, verbose=FALSE){

  # this is where we combine models and pad out eta 
  # with 0s as necessary to accommodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal

  # always collect if monitoring model is passed
  control$collect <- TRUE
  control$changes <- FALSE
  
  control$time.burnin <- 200
  control$time.interval <- 10
  control$time.samplesize <- 200
  
  # Must assign globals here, and after the optimization I should clean them up
  # and remove them from the environment.
  # The reason I must make them globals is because the cost function for
  # Bayesian Optimization can't have helper parameters (i think)
  global.control <<- control
  global.model <<- model
  global.model.mon <<- model.mon
  global.verbose <<- verbose
  global.current_ergm_state <<- ergm_state(nw, model=model.comb, proposal=proposal,
                                   stats=c(numeric(global.model$etamap$etalength), global.model.mon$nw.stats - global.model.mon$target.stats))
  global.current_best_dist <<- Inf
  
  path <- "/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptMore/"
  methods <- c("EI", "AEI", "UCB", "CL")
  seed <- 0
  times <- list(ei = numeric(), aei=numeric(), ucb=numeric(), cl=numeric())
  for(method in methods) {
    current_path <- paste0(path, method, "/")
    
    for (i in 1:2) {
      print("ITERATION")
      print(i)
      set.seed(seed)
      run <- NULL
      
      if(method == "EI") {
        start_time <- Sys.time()
        run <- runEI()
        times[["ei"]] <- c(times[["ei"]], Sys.time()-start_time)
      }
      else if (method == "AEI") {
        run <- runAEI()
        times[["aei"]] <- c(times[["aei"]], Sys.time()-start_time)
      }
      else if (method == "UCB") {
        run <- runUCB()
        times[["ucb"]] <- c(times[["ucb"]], Sys.time()-start_time)
      }
      else {
        run <- runCL()
        times[["cl"]] <- c(times[["cl"]], Sys.time()-start_time)
      }
      write.csv(as.data.frame(run$opt.path), file=paste0(current_path, method, "-", as.character(i), ".csv"), row.names=FALSE)
      
      seed <- seed + 1
    }
  }
  print(times)
  return()
  
  #run <- runPlain()
  #return()
  

  #print(run)
  #print(run$opt.path)
  #print(as.data.frame(run$opt.path))
  
  csv_files <- list.files("/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/CL/6x4/matern3_2+log/")
  max_number <- 0
  if (length(csv_files) != 0) {
    numbers <- as.numeric(gsub("CL-(\\d+)\\.csv", "\\1", csv_files))
    max_number <- max(numbers, na.rm = TRUE)
  }
  new_filename <- paste0("/Users/chris/Documents/uni_work/tergm_data/2d/2dBayesOptFixed/CL/6x4/matern3_2+log/CL-", max_number + 1, ".csv")
  write.csv(as.data.frame(run$opt.path), file=new_filename, row.names=FALSE)

  return()
  
  
  #write.csv(as.data.frame(run$opt.path), "opt_path_big_test.csv", row.names=FALSE)
  
  plot(run)
  #autoplot(obj.fun, render.levels = TRUE, show.optimum = TRUE) + geom_text(data = as.data.frame(res$opt.path), mapping = aes(label = dob), color = "white")
  
  print(run)
  return()
  newnetwork <- as.network(
                  ergm_state(global.nw, model=model.comb, proposal=proposal,
                             stats=c(numeric(global.model$etamap$etalength), 
                             global.model.mon$nw.stats - global.model.mon$target.stats))
                )
  
  theta <- unlist(run$x)
  eta <- ergm.eta(theta, global.model$etamap)
  
  print(eta)
  
  # remove all of the global from the environment
  #rm(global.control)
  #rm(global.model)
  #rm(global.model.mon)
  #rm(global.verbose)
    
  return(
    list(#nw.diff=z$nw.diff,
         newnetwork=newnetwork,
         eta=eta
         #history as well
    )
  )
  
}

runEI <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)
  
  parallelStop()
  
  mboControl = makeMBOControl()
  mboControl = setMBOControlTermination(mboControl, iters = 25)
  mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritEI())

  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  lrn = makeLearner("regr.km", predict.type = "se", nugget.estim = TRUE)  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  parallelStop()
  
  return(run)
}

runAEI <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)
  
  parallelStop()
  
  mboControl = makeMBOControl()
  mboControl = setMBOControlTermination(mboControl, iters = 25)
  mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritAEI())
  
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  lrn = makeLearner("regr.km", predict.type = "se", nugget.estim = TRUE)  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)

  return(run)
}

runUCB <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)
  parallelStop()
  
  mboControl = makeMBOControl()
  mboControl = setMBOControlTermination(mboControl, iters = 25)
  mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritCB())
  
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  lrn = makeLearner("regr.km", predict.type = "se", nugget.estim = TRUE)  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  
  return(run)
}

runCL <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)
  
  mboControl = makeMBOControl(propose.points = 4)
  mboControl = setMBOControlMultiPoint(mboControl, method = "cl", cl.lie=min)
  mboControl = setMBOControlTermination(mboControl, iters = 6)
  
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  lrn = makeLearner("regr.km", predict.type = "se", nugget.estim = TRUE)  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  parallelStop()
  
  return(run)
  
}


runPlain <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)
  
  #mboControl = makeMBOControl()
  #mboControl = setMBOControlTermination(mboControl, iters = 25)
  
  mboControl = makeMBOControl(propose.points = 4)
  mboControl = setMBOControlMultiPoint(mboControl, method = "cl", cl.lie=min)
  mboControl = setMBOControlTermination(mboControl, iters = 6)
  
  #mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritAEI(aei.use.nugget = TRUE))
  #mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritEI())
  #mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritCB())
  
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  lrn = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", nugget.estim = TRUE)  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  parallelStop()
  
  return(run)
  
}

runWithMultipleSamplesOfObjective <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, level="custom.objective")
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = parallelOptimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-6, upper=-4),
      makeNumericParam("theta2", lower=1, upper=3)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 3, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  
  des$y = apply(des, 1, obj.fun)
  
  mboControl = makeMBOControl()
  mboControl = setMBOControlTermination(mboControl, iters = 1)
  #mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritAEI(aei.use.nugget = TRUE))
  mboControl = setMBOControlInfill(mboControl, crit = makeMBOInfillCritEI())
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  
  parallelStop()
  
  return(run)
}

runWithMuleipleSamplesOfAcquisitionParallelDesign <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4, show.info = TRUE)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-10, upper=10),
      makeNumericParam("theta2", lower=-10, upper=10)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 25, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = parallelMap(optimCostFunctionWrapper, des[,1], des[,2], simplify = TRUE)

  
  mboControl = makeMBOControl(propose.points = 4)
  mboControl = setMBOControlMultiPoint(mboControl, method = "cl", cl.lie=min)

  mboControl = setMBOControlTermination(mboControl, iters = 6)
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  
  parallelStop()
  
  return(run)
}

runWithMultipleSamplesDesignMultipleAcquisition <- function() {
  parallelRegisterLevels(levels = "objective")
  parallelStart("multicore", 4)
  
  obj.fun <- makeSingleObjectiveFunction(
    name = "optimCostFunction",
    fn = optimCostFunction,
    par.set = makeParamSet(
      makeNumericParam("theta1", lower=-6, upper=-4),
      makeNumericParam("theta2", lower=1, upper=3)
    ),
    minimize = TRUE,
    noisy = TRUE
  )
  
  des = generateDesign(n = 20, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  des$y = apply(des, 1, parallelOptimCostFunction)
  
  mboControl = makeMBOControl(propose.points = 4)
  mboControl = setMBOControlMultiPoint(mboControl, method = "cl", cl.lie = min)
  mboControl = setMBOControlTermination(mboControl, iters = 20)
  lrn = makeMBOLearner(mboControl, obj.fun, nugget.estim = TRUE)
  
  run = mbo(obj.fun, design = des, learner = lrn, control = mboControl, show.info = TRUE)
  
  parallelStop()
  
  return(run)
}


mahalanobisDist1 <- function(target_statistics) {
  
  cov_matrix <- cov(target_statistics)
  cov_out <<- cov_matrix

  # check for singular matrix
  if (rcond(cov_matrix) < .Machine$double.eps) {
    return(10^10)
  }

  inv_cov_matrix <- solve(cov_matrix)
  
  
  #dist <- sqrt(mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix))
  #dist <- log(mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix))
  dist <- mahalanobis(colMeans(target_statistics), rep(0, ncol(target_statistics)), inv_cov_matrix, inverted=TRUE)
  if (dist < 0) {
    return(10^10)
  }
  return(dist)
}

optimCostFunction <- function(theta) {
  #theta <-c(-5.194653, 2.194362)
  #theta1 <<- theta
  
  # call many times
  # re-run every time update theta
  eta <- ergm.eta(theta, global.model$etamap)
  eta.comb <- c(eta, numeric(global.model.mon$etamap$etalength))

  z <- tergm_MCMC_slave(global.current_ergm_state, eta.comb, global.control, global.verbose)
  
  target_statistics <- z$statsmatrix[,-seq_len(nparam(global.model, canonical = TRUE))]
  #plot(density(target_statistics[,1]))
  
  #print(target_statistics)
  # optional to change ergm state every time
  # global.current_ergm_state <<- z$state
  
  res <- mahalanobisDist1(target_statistics)
  
  # if (global.current_best_dist > res) {
  #   global.current_best_dist <<- res
  #   global.current_ergm_state <<- z$state
  # }
  

  #targets <<- target_statistics
  
  # For debugging
  # mahalanob <<- mahalanobisDist(target_statistics)
  #return(log(1+res))
  #return(sqrt(1+res))
  return(res)
  #return(mahalanobisDist(target_statistics))  
}

optimCostFunctionWrapper <- function(...) {
  theta <- c(...)
  return(optimCostFunction(theta))
}

parallelOptimCostFunction <- function(theta) {
  rep.theta = replicate(4, theta, simplify = FALSE)
  
  result = parallelMap(optimCostFunction, rep.theta, simplify = TRUE)

  return(mean(result))
}






