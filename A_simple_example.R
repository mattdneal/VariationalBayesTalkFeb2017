# A simple example

iterate <- function(par, hyperpar, X) {
  l0 <- hyperpar$l0
  m0 <- hyperpar$m0
  a0 <- hyperpar$a0
  b0 <- hyperpar$b0
  N <- length(X)
  
  out <- list()
  
  # constant across iterations
  out$mN <- (l0*m0 + sum(X)) / (l0 + N)
  out$aN <- a0 + (N + 1) / 2
  
  # First update bN using the last value of lN
  out$bN <- b0 + 1 / 2 * ((l0 + N) * (1 / par$lN + out$mN^2)
                          - 2 * (l0 * m0 + sum(X) * out$mN)
                          + sum(X^2)
                          + l0 * m0^2)
  
  # Now update lN with the new value for bN
  out$lN <- (l0 + N) * out$aN / out$bN
  return(out)
}

q_mean <- function(x, par) {
  dnorm(x, mean=par$mN, sd=1 / sqrt(par$lN))
}

q_precision <- function(x, par) {
  dgamma(x, shape = par$aN, rate=par$bN)
}

q <- function(x, par) {
  as.numeric(apply(x, 1, function(row) q_mean(row[1], par) * q_precision(row[2], par)))
}

plot_posterior <- function(par, main="Posterior") {
  
  nrange <- qnorm(c(0.005, 0.995), mean=par$mN, sd=1 / sqrt(par$aN/par$bN * par$lN))
  nrange <- seq(nrange[1], nrange[2], length.out=1000)
  
  grange <- qgamma(c(0.005, 0.995), shape = par$aN, rate=par$bN)
  grange <- seq(grange[1], grange[2], length.out=1000)
  
  layout(t(1:2))
  
  plot(nrange, q_mean(nrange, par), type="l", xlab="", ylab="", main=paste(main, "mean distribution"))
  plot(grange, q_precision(grange, par), type="l", xlab="", ylab="", main=paste(main, "precision distribution"))
  
  layout(1)
}

plot_contour <- function(par) {
  
  nrange <- qnorm(c(0.005, 0.995), mean=par$mN, sd=1 / sqrt(par$aN/par$bN * par$lN))
  nrange <- seq(nrange[1], nrange[2], length.out=100)
  
  grange <- qgamma(c(0.005, 0.995), shape = par$aN, rate=par$bN)
  grange <- seq(grange[1], grange[2], length.out=100)
  coords <- as.matrix(expand.grid(nrange, grange))
  z <- matrix(q(coords, par), 100)
  contour(x=nrange, y=grange, z=z, xlab="Mean", ylab="Precision")
}

plot_prior <- function(hyperpar) {
  par <- list()
  par$mN <- hyperpar$m0
  par$aN <- hyperpar$a0
  par$bN <- hyperpar$b0
  par$lN <- hyperpar$l0
  plot_posterior(par, main="Prior")
}

X <- rnorm(10, mean = 15, sd=5)
hyperpar <- list()
hyperpar$l0 <- 0.001
hyperpar$m0 <- 0
hyperpar$a0 <- 1
hyperpar$b0 <- 1

plot_prior(hyperpar)

par <- list()
# These will be calculated in the first iteration
par$mN <- NA
par$aN <- NA
par$bN <- NA
# But we need a starting point for lN
par$lN <- 1

for (i in 1:10) {
  par <- iterate(par, hyperpar, X); print(par)
  plot_contour(par)
  readline()
}
plot_contour(par)
  


X <- rnorm(1000, mean = 15, sd=5)
hyperpar <- list()
hyperpar$l0 <- 0.001
hyperpar$m0 <- 0
hyperpar$a0 <- 1
hyperpar$b0 <- 1

plot_prior(hyperpar)

par <- list()
# These will be calculated in the first iteration
par$mN <- NA
par$aN <- NA
par$bN <- NA
# But we need a starting point for lN
par$lN <- 1

for (i in 1:10) {
  par <- iterate(par, hyperpar, X); print(par)
  plot_contour(par)
  readline()
}
plot_contour(par)
