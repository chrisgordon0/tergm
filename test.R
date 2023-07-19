library(tsna)
library(ndtv)
library(htmlwidgets)
library(latticeExtra)
#library(tergm)
library(mlrMBO)

library(ergm)
library(ergm.multi)

library(devtools)
library(ggplot2)

library(foreach)
library(doParallel)
library(parallelMap)
library(GPFDA)

load_all()


simple_test <- function() {
  n<-20
  do.plot <- FALSE
  g0<-network.initialize(n,dir=FALSE)
  
  #                    edges, mean.age
  target.stats<-c(     n*1/2,       10)
  
  coef.exact<-function(density,duration)
    list(form=-log(((1+exp(logit(1-1/duration)))/(density/(1-density)))-1),
         diss=logit(1-1/duration))
  
  
  truth <- coef.exact(target.stats[1]/network.dyadcount(g0),
                      target.stats[2])
  
  # Get a deliberately bad starting network.
  set.seed(3)
  g1<-san(g0~meandeg,target.stats=target.stats[1],verbose=TRUE)
  
  # Fit the model with very poor starting values.
  set.seed(4)
  #dynfit<-tergm(g1 ~ Form(~edges), targets=~edges+mean.age, estimate="EGMME",target.stats=target.stats[-3], constraints=~., verbose=FALSE, control=control.tergm(SA.plot.progress=do.plot,SA.phase2.levels.min=2, SA.phase2.levels.max=4, SA.phase2.repeats = 10, SA.restart.on.err=FALSE,init=c(-log(.95/.05))))
  dynfit<-tergm(g1 ~ Form(~edges) + Persist(~edges), targets=~edges+mean.age, estimate="EGMME",target.stats=target.stats[-3], constraints=~., verbose=FALSE,control=control.tergm(SA.plot.progress=do.plot,SA.phase2.levels.min=2, SA.phase2.levels.max=4, SA.phase2.repeats = 10, SA.restart.on.err=FALSE,init=c(-log(.95/.05), 1)))
  #print(dynfit)
  #return(dynfit[1])
  
  #expect_warning(expect_error(print(summary(dynfit)), NA), NA)
  #expect_warning(expect_error(mcmc.diagnostics(dynfit), NA), NA)
  
  #expect_equal(unlist(truth),coef(dynfit),tolerance=0.05,ignore_attr=TRUE)
}
start_time <- Sys.time()
simple_test()
end_time <- Sys.time()

print("TOTAL EXECUTION TIME:")
print(end_time - start_time)

# I think I want from -6 to 0 and from 0 to 6 


plot_objective <- function () {
  # Create a grid of theta1 and theta2 values
  theta1 <- seq(-6, 0, 0.1)
  theta2 <- seq(0, 6, 0.1)

  num_theta1 <- length(theta1)
  num_theta2 <- length(theta2)

  
  # Create an empty matrix to store the objective function values
  objective_matrix <- matrix(NA, nrow = num_theta1, ncol = num_theta2)
  
  
  # Loop over the theta1 and theta2 values and calculate the objective function values
  for (i in 1:num_theta1) {
    for (j in 1:num_theta2) {
      
      theta_value <<- c(theta1[i], theta2[j])
      
      result <- simple_test()
      result <- unlist(result)

      objective_matrix[i, j] <- result
    }
  }

  write.csv(objective_matrix, file = "objective_matrix.csv")
  
  dataframe_data <- as.data.frame(objective_matrix)
  
  # Save the dataframe to a CSV file
  write.csv(dataframe_data, file = "objective_dataframe.csv", row.names = FALSE)
}

plotObjectiveMatrix <- function(file_path) {
  library('plot.matrix')
  # library('plot3D')
  
  # Load the objective matrix from the CSV file
  objective_df <- read.csv(file_path, header = TRUE)
  objective_mat <- data.matrix(objective_df, rownames.force = NA)
  contour(x = seq(-6, 0, 0.1), y = seq(0, 6, 0.1), objective_mat, levels=c(1,2,3,4,5,10,15,20,30,40,50,100,1000,200,3000))
  # persp3D(z = objective_mat, theta = 120)
  
}

plotObjectiveMatrixSquared <- function(file_path) {
  library('plot.matrix')
  # library('plot3D')
  
  # Load the objective matrix from the CSV file
  objective_df <- read.csv(file_path, header = TRUE)
  objective_mat <- data.matrix(objective_df, rownames.force = NA)
  
  fun <- function(x) {
    x^2
  }
  
  objective_mat <- apply(objective_mat, c(1, 2), fun)
  
  contour(x = seq(-6, 0, 0.1), y = seq(0, 6, 0.1), objective_mat, levels=c(1,2,3,4,5,10,15,20,30,40,50,100,1000,200,3000,5000))
  # persp3D(z = objective_mat, theta = 120)
}


plotDensityOf100Samples <- function() {
  data <- c(7.719301, 2.284666, 3.494724, 6.617824, 1.350646, 6.734550, 1.447522, 5.791716, 3.415927, 8.695420,
            3.567117, 4.139120, 1.266246, 2.663978, 2.049393, 3.533716, 1.385487, 5.019851, 3.282313, 2.560116,
            4.620162, 2.886019, 3.066905, 1.660840, 5.205677, 3.093498, 3.786891, 2.404962, 1.603466, 1.675048,
            1.548542, 3.570324, 2.010171, 2.236274, 5.223538, 3.050817, 5.914710, 2.636048, 2.377787, 4.360276,
            2.289258, 3.172092, 1.978822, 2.489084, 4.064381, 3.115198, 10.138987, 1.624544, 2.563657, 3.175003,
            5.769803, 9.889008, 5.277521, 5.971271, 6.975253, 1.642959, 1.921086, 4.187314, 3.032391, 6.518522,
            3.070071, 5.926540, 2.287304, 2.948136, 3.690434, 3.601427, 6.265596, 1.564738, 3.790398, 2.444413,
            1.450112, 2.662356, 1.601634, 3.421766, 3.554014, 1.896244, 1.257736, 9.318663, 5.832165, 6.923757,
            2.764944, 1.636017, 5.826695, 3.169947, 4.239614, 5.868942, 3.944731, 4.086817, 5.430511, 3.713535,
            1.932094, 2.811196, 7.890121, 9.391567, 2.590621, 1.532060, 1.738207, 2.229768, 1.461112, 2.107796)
  
  plot(data, type = "o", pch = 16, col = "blue", xlab = "Index", ylab = "Value")
  
  density_est <- density(data)
  plot(density_est, main = "Density Plot of Data", xlab = "Value", ylab = "Density")
}

plotOptPath <- function() {
  optPath <- read.csv("opt_path_big_test.csv", header = TRUE)
  x <- seq(1, length(optPath$y), 1)
  print(optPath$y)
  plot(x, log(optPath$y), type="l", ylab="log(Objective Function)", xlab="iterations")
  abline(v=40, col='red')
  lines(lowess(x, log(optPath$y)))
}
#plotOptPath()
# start_time <- Sys.time()
# simple_test()
# end_time <- Sys.time()
# 
# print("TOTAL TIME:")
# print(end_time - start_time)


#plot_objective()
#plotObjectiveMatrix("objective_dataframe.csv")
#plotObjectiveMatrixSquared("objective_dataframe.csv")
#plotDensityOf100Samples()

# n = 2
# values = rep(0,n)
# for (i in 1:n) {
#   result <- simple_test()
#   result <- unlist(result)
# 
#   values[i] = result
# }
# 
# print(values)
# print(mean(values))
# print(var(values))



# data(samplk)
# 
# samp.series <- NetSeries(list(samplk1,samplk2,samplk3))
# 
# dat <- data(florentine)
# summary(dat)
# 
# theta.persist <- log(9)
# 
# data(florentine)
# 
# set.seed(1)
# egmme.fit <- tergm(
#   flobusiness ~
#     Form(~ edges + gwesp(0, fixed=T)) +
#     Persist(~ offset(edges)),
#   targets = ~ edges + gwesp(0, fixed=T),
#   offset.coef = theta.persist,
#   estimate = "EGMME",
#   control = control.tergm(SA.plot.progress=TRUE)
# )
# 
# egmme.fit
# 
# mcmc.diagnostics(egmme.fit, which="plots") # only returns the plots
