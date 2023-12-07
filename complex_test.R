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
library(threejs)
library(statnet.common)

load_all()

library(ergm)
library(tergm)

n.nodes = 1000
mean.deg = 2

#duration = seq(5, 50, by = 5)
duration = c(5)

#nsim = 10
nsim = 1

for(n in 1:nsim) {
  for(i in length(duration):1){
    dur = duration[i]
    gamma = log(dur -1)
    dyn.fit <- try(tergm(g ~ Form(~edges + degree(1:2)+ gwesp(0.5, fixed = T)) + Diss(~offset(edges)),
                          estimate = "EGMME", 
                          targets = ~edges+degree(1:2)+gwesp(0.5, fixed = T),
                          offset.coef = gamma, control=control.tergm(EGMME.MCMC.burnin.max = 10000)))#, control = control.stergm(init.form = theta0)
    print(dyn.fit)
  }
  print(paste("Completing simulation", n, "of", nsim))
}

#save(coef.diff, mean.diff, sd.diff, xsec.coef, meanstats.dev, file = "1000node_Tri_stoch.RData")
