library(GPFDA)
input <- matrix(c(1,2,3,4,5,1,2,3,4,5), ncol=5, nrow=2, byrow=TRUE)
response <- c(1,2)

fit <- gpr(input=input, response=response)
fit
