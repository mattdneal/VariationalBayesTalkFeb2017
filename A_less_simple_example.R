#install.packages("MASS")

generate_gmm_data <- function(means, pi, n) {
  if (nrow(means) != length(pi)) stop("Invalid parameters")
  Z <- matrix(0, n, length(pi)) 
  X <- matrix(0, n, ncol(means))
  dist <- as.numeric(rmultinom(1, n, pi))
  row <- 0
  lambda <- list()
  for (k in 1:length(pi)) {
    lambda[[k]] <- rWishart(1, ncol(means), diag(1, ncol(means)))[,, 1]
    if (dist[k] > 0) {
      for (i in 1:dist[k]) {
        row <- row + 1
        Z[row, k] <- 1
        X[row, ] <- MASS::mvrnorm(1, means[k, ], lambda[[k]])
      }
    }
  }
  return(list(X=X, Z=Z, lambda=lambda, means=means))
}

plot_gmm <- function(gmm) {
  col <- apply(gmm$Z, 1, which.max)
  plot(gmm$X, col=col)
}

iterate <- function(par, hyperpar, X) {
  K <- hyperpar$K
  D <- ncol(X)
  N <- nrow(X)
  v <- par$v
  out <- list()
  
  W0.inv <- solve(hyperpar$W0)
  #first we calculate r
  E_lambda <- numeric(K)
  E_pi <- numeric(K)
  
  for (k in 1:K) {
    for (i in 1:D) {
      E_lambda[k] <- E_lambda[k] + digamma((v[k] + 1 - i) / 2)
    }
    E_lambda[k] <- E_lambda[k] + D * log(2) + log(det(par$W[[k]]))
    E_pi[k] <- digamma(par$a[k]) - digamma(sum(par$a))
  }
  
  E_lambda <- exp(E_lambda)
  E_pi <- exp(E_pi)
  
  r <- matrix(0, N, K)
  for (n in 1:N) {
    for (k in 1:K) {
      r[n, k] <- E_pi[k] * sqrt(E_lambda[k]) * exp(-D / (2 * par$b[k]) - par$v[k] / 2 * t(X[n, ] - par$m[k, ]) %*% par$W[[k]] %*% (X[n, ] - par$m[k, ]))
    }
  }
  
  plot(X, col=apply(r, 1, which.max), main="r based on input par")
  points(par$m, pch=2, col=1:nrow(par$m))
  for (i in 1:nrow(par$m)) {
    points(MASS::mvrnorm(1000, par$m[i,], solve(par$W[[i]] * par$v[i])), pch=3, col=i)
  }
  r_rowsums <- rowSums(r)
  
  for (n in 1:N) {
    r[n, ] <- r[n, ] / r_rowsums[n]
  }
  out$r <- r
  #Now we need to update the other parameters
  
  r_colsums <- colSums(r) # N_k
  x.bar <- matrix(0, K, D)
  for (k in 1:K) {
    x.bar[k, ] <- 1 / r_colsums[k] * colSums(r[, k] * X)
  }
  S <- list()
  for (k in 1:K) {
    S[[k]] <- matrix(0, D, D)
    for (n in 1:N) {
      S[[k]] <- S[[k]] + r[n, k] * (X[n, ] - x.bar[k, ]) %*% t(X[n, ] - x.bar[k, ])
    }
    S[[k]] <- 1 / r_colsums[k] * S[[k]]
  }
  
  b <- hyperpar$b0 + r_colsums
  m <- 1 / b * (hyperpar$b0 * hyperpar$m0 + r_colsums * x.bar )
  W <- list()
  for (k in 1:K) {
    W.inv.temp <- (W0.inv 
                   + r_colsums[k] * S[[k]] 
                   + hyperpar$b0 * r_colsums[k] 
                   / (hyperpar$b0 + r_colsums[k]) 
                   * (x.bar[k, ] - hyperpar$m0) %*% t(x.bar[k, ] - hyperpar$m0))
    W[[k]] <- chol2inv(chol(W.inv.temp))
  }
  a <- hyperpar$a0 + r_colsums
  v <- hyperpar$v0 + r_colsums
  
  out$a <- a
  out$m <- m
  out$W <- W
  out$b <- b
  out$v <- v
  plot(X, col=apply(r, 1, which.max), main="output")
  points(out$m, pch=2, col=1:nrow(out$m))
  for (i in 1:nrow(out$m)) {
    points(MASS::mvrnorm(1000, out$m[i,], solve(out$W[[i]] * out$v[i])), pch=3, col=i)
  }
  return(out)
}

plot_iteration <- function(X, it) {
  col <- apply(it$r, 1, which.max)
  plot(X, col=col)
}

out <- generate_gmm_data(matrix(rnorm(8)*3, 4), 1:4 / 10, 1000)
#out <- generate_gmm_data(matrix(rnorm(4)*2, 2), 1:2 / 3, 1000)
plot_gmm(out)
X <- out$X
K <- ncol(out$Z)

#set some prior hyperparameters
hyperpar <- list()
hyperpar$K <- K
hyperpar$a0 <- 0.1
hyperpar$v0 <- 2
hyperpar$b0 <- 0.01
hyperpar$m0 <- c(0,0)
hyperpar$W0 <- diag(1, 2)

km <- kmeans(X, hyperpar$K)

#make a first guess at our parameters
par <- list()
par$a <- rep(hyperpar$a0, hyperpar$K)
par$v <- rep(hyperpar$v0, hyperpar$K)
par$m <- matrix(rnorm(length(out$means)), hyperpar$K) #km$centers
par$W <- list()
for (k in 1:hyperpar$K) par$W[[k]] <- hyperpar$W0
par$b <- rep(1, hyperpar$K)

for (i in 1:100) {
  par <- iterate(par, hyperpar, X)
  plot_iteration(X, par)
  print(par$m)
  readline()
}
