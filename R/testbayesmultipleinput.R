library(mlrMBO)

f <- function() {

  obj.fun <- makeSingleObjectiveFunction(
    name = "my_sphere",
    fn = function(x) {
      -sum(x*x)
    },
    par.set = makeParamSet(
      makeNumericVectorParam("x", len = 2L, lower = -5, upper = 5)
    ),
    minimize = TRUE
  )
  
  des = generateDesign(n = 5, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
  
  des$y = apply(des, 1, obj.fun)
  
  surr.km = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", control = list(trace = FALSE))
  
  control = makeMBOControl()
  control = setMBOControlTermination(control, iters = 10)
  control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())
  
  
  run = mbo(obj.fun, design = des, learner = surr.km, control = control, show.info = TRUE)
  
  
  print(run)

}