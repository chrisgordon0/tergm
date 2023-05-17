library(tsna)
library(ndtv)
library(htmlwidgets)
library(latticeExtra)
#library(tergm)
library(rBayesianOptimization)

library(ergm)
library(ergm.multi)

library(devtools)

load_all()

# lookup devtools - devmode
# debug function in R

# Example 1: Optimization
## Set Pred = 0, as placeholder

#Test_Fun <- function(x) {
#  list(Score = -x^2, Pred = 0)
#}

## Set larger init_points and n_iter for better optimization result
#OPT_Res <- BayesianOptimization(Test_Fun,
#                                bounds = list(x = c(0, 1)),
#                                init_points = 2, n_iter = 2,
#                                acq = "ucb", kappa = 1, eps = 0.0,
#                                verbose = TRUE)

#OPT_Res$Best_Par
#OPT_Res$Best_Value

set.seed(1)


data(samplk)

samp.series <- NetSeries(list(samplk1,samplk2,samplk3))

dat <- data(florentine)
summary(dat)

theta.persist <- log(9)

data(florentine)

set.seed(1)
egmme.fit <- tergm(
  flobusiness ~ 
    Form(~ edges + gwesp(0, fixed=T)) + 
    Persist(~ offset(edges)),
  targets = ~ edges + gwesp(0, fixed=T),
  offset.coef = theta.persist,
  estimate = "EGMME",
  control = control.tergm(SA.plot.progress=TRUE)
)

egmme.fit

mcmc.diagnostics(egmme.fit, which="plots") # only returns the plots
