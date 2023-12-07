#library(tergm)
library(coda)
library(parallel)
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

nw <- network.initialize(1000, directed = FALSE)

edges_target <- 1000
degree1_target <- 200
degree2_target <- 350
gwesp_targets <- 3*c(1, seq(from = 10, to = 100, by = 10))

# target durations
mean_duration_targets <- c(15, 50, 100)

simulation_style <- c("old", "new")

ncores <- 40
nsteps_burnin <- 4000
nsteps_sample <- 36000

targets_grid = expand.grid(edges = edges_target,
                           degree1 = degree1_target,
                           degree2 = degree2_target,
                           gwesp = gwesp_targets,
                           mean_duration = mean_duration_targets)

index <- 1

start_time <- Sys.time()

set.seed(2)

for (i in 1:10) {
  set.seed(i)
  terg <- tergm(nw ~ Form(~edges + degree(1) + degree(2) + gwesp(0.5, fixed = TRUE)) + Persist(~edges),
                targets = ~edges + degree(1) + degree(2) + gwesp(0.5, fixed = T) + mean.age,
                target.stats = unname(unlist(targets_grid[index,])),
                estimate = "EGMME") 
}

end_time <- Sys.time()
print(end_time-start_time)
print(coef(terg))
print(end_time-start_time)


estimate_model <- function(index) {
  library(ergm)

  ergm_fit <- ergm(nw ~ edges + degree(1) + degree(2) + gwesp(0.5, fixed = TRUE),
                   target.stats = unname(unlist(targets_grid[index,])),
                   control = list(init.MPLE.samplesize = network.dyadcount(nw),
                                  MCMLE.effectiveSize = NULL,
                                  MCMC.samplesize = 5e3,
                                  MCMC.burnin = 5e3,
                                  MCMC.interval = 5e3,
                                  SAN.nsteps = 5e7))

  coef(ergm_fit)
}

# estimate models in parallel
cl <- makeCluster(ncores)

clusterExport(cl, "nw")
clusterExport(cl, "targets_grid")

st <- Sys.time()
rv <- parLapply(cl, seq_len(NROW(targets_grid)), estimate_model)
et <- Sys.time()
print(et - st)
# Time difference of 1.611947 mins

stopCluster(cl)

ergm_coefs <- do.call(rbind, rv)

argsgrid <- expand.grid(edges = edges_target,
                        degree1 = degree1_target,
                        degree2 = degree2_target,
                        gwesp = gwesp_targets,
                        mean_duration = mean_duration_targets,
                        style = simulation_style)

argsgrid <- cbind(argsgrid, ergm_coefs)

simulate_model <- function(index) {
  library(tergm)
  library(coda)

  style <- argsgrid[index, 6]
  ergm_coef <- unlist(argsgrid[index, 7:10])
  mean_duration <- argsgrid[index, 5]

  sim_formula <- nw ~ Form(~edges + degree(1) + degree(2) + gwesp(0.5, fixed = TRUE)) + Persist(~edges) # use this formula
  if(style == "new") {
    sim_coef <- c(ergm_coef[1] - log(mean_duration), ergm_coef[-1], log(mean_duration - 1))
  } else {
    sim_coef <- c(ergm_coef[1] - log(mean_duration - 1), ergm_coef[-1], log(mean_duration - 1))
  }

  rv <- simulate(sim_formula,
                 coef = sim_coef,
                 time.burnin = nsteps_burnin,
                 time.slices = nsteps_sample,
                 control = list(MCMC.burnin.min = 1e5, MCMC.burnin.max = 1e7),
                 output = "stats",
                 monitor = ~edges + degree(1) + degree(2) + gwesp(0.5, fixed = TRUE),
                 dynamic = TRUE)

  m <- colMeans(rv)
  ess <- effectiveSize(rv)
  se <- apply(rv, 2, sd)/sqrt(ess)

  list(m = m,
       ess = ess,
       se = se,
       index = index)
}

# run sims in parallel
cl <- makeCluster(ncores)

clusterExport(cl, "nw")
clusterExport(cl, "argsgrid")
clusterExport(cl, "nsteps_burnin")
clusterExport(cl, "nsteps_sample")

st <- Sys.time()
rv <- parLapply(cl, seq_len(NROW(argsgrid)), simulate_model)
et <- Sys.time()
print(et - st)
# Time difference of 13.42398 hours

stopCluster(cl)

## include sim means in output object
ag <- argsgrid

m <- do.call(rbind, lapply(rv, `[[`, "m"))
colnames(m) <- paste0("sim_mean_", c("edges", "degree1", "degree2", "gwesp"))
ess <- do.call(rbind, lapply(rv, `[[`, "ess"))
colnames(ess) <- paste0("sim_ess_", c("edges", "degree1", "degree2", "gwesp"))
se <- do.call(rbind, lapply(rv, `[[`, "se"))
colnames(se) <- paste0("sim_se_", c("edges", "degree1", "degree2", "gwesp"))
ag <- cbind(ag, m, ess, se)

save(ag, file = "gwesp_data.rdata")
