library(GPFDA)
require(MASS)

set.seed(123)
nrep <- 30
n1 <- 25
n2 <- 25
n3 <- 25
N <- 3
n <- n1+n2+n3
input1 <- sapply(1:n1, function(x) (x - min(1:n1))/max(1:n1 - min(1:n1)))
input2 <- input1
input3 <- input1

# storing input vectors in a list
Data <- list()
Data$input <- list(input1, input2, input3)

# true hyperparameter values
nu0s <- c(6, 4, 2)
nu1s <- c(0.1, 0.05, 0.01)
a0s <- c(500, 500, 500)
a1s <- c(100, 100, 100)
sigm <- 0.05
hp <- c(nu0s, log(nu1s), log(a0s), log(a1s), log(sigm))

# Calculate covariance matrix
Psi <- mgpCovMat(Data=Data, hp=hp)

ns <- sapply(Data$input, length)
idx <- c(unlist(sapply(1:N, function(i) rep(i, ns[i])))) 


mu <- c( 5*input1, 10*input2, -3*input3)
Y <- t(mvrnorm(n=nrep, mu=mu, Sigma=Psi))
response <- list()
for(j in 1:N){
  response[[j]] <- Y[idx==j,,drop=F]
}
# storing the response in the list
Data$response <- response


res <- mgpr(Data=Data, m=20, meanModel = 't')
