rm(list = ls())

library(coda)
library(rstan)
library(EBRR)

set.seed(372981)

dir <- ""

dataname <- "prostate" # options: "prostate" or "glucose"

burn <- 1000
num.samp <- 1000
chains <- 10
cores <- 2


if (dataname == "prostate") {
  # Got from here: https://cran.r-project.org/src/contrib/Archive/lasso2/
  load(paste(dir, "Data/Prostate.rda", sep = ""))

  y <- Prostate[, "lpsa"]
  y <- y
  X <- as.matrix(Prostate[, names(Prostate) != "lpsa"])
  n2 <- ncol(X)
  m <- nrow(X)
  X <- apply(X, 2, function(x) {(x - mean(x))/sd(x)}) # Original paper standardizes
  y <- (y - mean(y))/sd(y)
  x <- X
} else if (dataname == "glucose") {
  data <- read.csv(paste(dir, "Data/Shahetal_Metab_Dataset.csv", sep = ""),
                   stringsAsFactors = FALSE)
  data <- apply(data, c(1, 2), function(x) ifelse(x == ".", NA, x))
  data[, "sex"] <- ifelse(data[, "sex"] == "M", 1, 0)
  data <- data.frame(cbind(apply(data, 2, function(x) {as.numeric(x)})))
  data$htn <- ifelse(data$htn == 2, 0, 1)
  data$dyslipid <- ifelse(data$dyslipid == 2, 0, 1)
  data$cad <- ifelse(data$cad == 2, 0, 1)
  data$dm <- ifelse(data$dm == 2, 0, 1)
  # Construct family variables, should regress out
  data$family2 <- ifelse(data$family == 2, 1, 0)
  data$family3 <- ifelse(data$family == 3, 1, 0)
  data$family4 <- ifelse(data$family == 4, 1, 0)
  data$family5 <- ifelse(data$family == 5, 1, 0)
  data$family6 <- ifelse(data$family == 6, 1, 0)
  data$family7 <- ifelse(data$family == 7, 1, 0)
  data$family8 <- ifelse(data$family == 8, 1, 0)
  data$family <- NULL

  data <- data[complete.cases(data), ]
  data <- data[, apply(data, 2, var) != 0] # Kick out any variables that are constant
  # Either parametrization produces similar results
  y <- data$GLU  # Looks normal enough
  X <- as.matrix(data[, !colnames(data) %in% c("GLU", "Obs")])
  m <- length(y)
  Z <- cbind(rep(1, m), X[, grepl("fam", colnames(X))])
  P <- diag(m) - Z%*%solve(crossprod(Z))%*%t(Z)
  y <- P%*%y
  X <- P%*%X[, !grepl("fam", colnames(X))]
  # Not sure if I should standardize at first or later - I think later?
  y <- y/sd(y)
  X <- apply(X, 2, function(x) {(x)/sd(x)})
  n2 <- ncol(X)
  x <- X
  y <- c(y)
}

ebrr <- estRegPars(y = y, X = x, comp.q = FALSE, mom = FALSE)
sigma <- sqrt(ebrr$sigma.epsi.sq.hat)
tau <- sqrt(ebrr$sigma.beta.sq.hat)

stanmod1 <- stan_model(file =  paste(dir, "/Code/powerregfixlamq_noncentered.stan", sep = ""))
stanmod2 <- stan_model(file =  paste(dir, "/Code/powerregfixlamq_centered.stan", sep = ""))
stanmod3 <- stan_model(file =  paste(dir, "/Code/powerregfixlamq_naive.stan", sep = ""))

z2.start <- solve(crossprod(X)/sigma^2 + diag(n2)/tau^2)%*%crossprod(X, y)/sigma^2
vz2.start <- sigma^2*solve(crossprod(X)/sigma^2 + diag(n2)/tau^2)%*%crossprod(X)%*%solve(crossprod(X)/sigma^2 + diag(n2)/tau^2) + diag(n2)*10^(-18)

inits <- lapply(1:chains, function(chain) {
  list("z2" = c(z2.start + t(chol(vz2.start))%*%rnorm(n2)))
})

# For multiple q, try all three model specifications
qs <- seq(0.2, 1.8, by = 0.2)
for (i in 1:length(qs)) {
  q <- qs[i]
  lambda <- ((gamma(3/q)/gamma(1/q))/tau^2)^(q/2)

  # Make sure all starting values are the same
  inits1 <- lapply(inits, function(init, lambda, q) {
    z2 <- init[["z2"]]
    n2 <- length(z2)
    delta <- pi/2
    div <- ((2.0^(-0.5)*lambda^(-1.0/q)/sin(q*delta/2.0)^(1.0/2.0)/sin((2.0 - q)*delta/2.0)^((2.0 - q)/(2.0*q))*sin(delta)^(1.0/q)))
    list("lxi" = rep(0, n2), "logitdelta" = rep(0, n2),
         "w" = z2/div)
  }, lambda = lambda, q = q)

  inits2 <- lapply(inits, function(init) {
    z2 <- init[["z2"]]
    n2 <- length(z2)
    list("z2" = z2, "lxi" = rep(0, n2), "logitdelta" = rep(0, n2))
  })

  inits3 <- lapply(inits2, function(x, q, lambda) {list("w" = x[["z2"]]/lambda^(-1.0/q))},
                   q = q, lambda = lambda)

  fitf1 <- sampling(stanmod1,
                    data = list(m = m, n2 = n2, y = y, x = x,
                                lambda = lambda, q = q, sigma = sigma),    # named list of data
                    chains = chains,             # number of Markov chains
                    warmup = burn,          # number of warmup iterations per chain
                    iter = burn + num.samp,    # total number of iterations per chain
                    cores = cores,              # number of cores (could use one per chain)
                    refresh = 1,
                    init = inits1)
  res1 <- as.array(fitf1)

  fitf2 <- sampling(stanmod2,
                    data = list(m = m, n2 = n2, y = y, x = x,
                                lambda = lambda, q = q, sigma = sigma),    # named list of data
                    chains = chains,             # number of Markov chains
                    warmup = burn,          # number of warmup iterations per chain
                    iter = burn + num.samp,    # total number of iterations per chain
                    cores = cores,              # number of cores (could use one per chain)
                    refresh = 1, init = inits2)
  res2 <- as.array(fitf2)

  fitf3 <- sampling(stanmod3,
                    data = list(m = m, n2 = n2, y = y, x = x,
                                lambda = lambda, q = q, sigma = sigma),    # named list of data
                    chains = chains,             # number of Markov chains
                    warmup = burn,          # number of warmup iterations per chain
                    iter = burn + num.samp,    # total number of iterations per chain
                    cores = cores,              # number of cores (could use one per chain)
                    refresh = 1, init = inits3)
  res3 <- as.array(fitf3)

  v1s <- c(paste("lxi[", 1:n2, "]", sep = ""),
           paste("logitdelta[", 1:n2, "]", sep = ""),
           paste("w[", 1:n2, "]", sep = ""))
  v2s <- c(paste("lxi[", 1:n2, "]", sep = ""),
           paste("logitdelta[", 1:n2, "]", sep = ""),
           paste("z2[", 1:n2, "]", sep = ""))
  v3s <- c(paste("z2[", 1:n2, "]", sep = ""))

  ess1 <- apply(apply(res1[, , v1s], c(2, 3), function(x) {
    effectiveSize(x)
  }), 1, min)
  ess2 <- apply(apply(res2[, , v2s], c(2, 3), function(x) {
    effectiveSize(x)
  }), 1, min)
  ess3 <- apply(apply(res3[, , v3s], c(2, 3), function(x) {
    effectiveSize(x)
  }), 1, min)

  z21s <- res1[, , v3s]
  z22s <- res2[, , v3s]
  z23s <- res3[, , v3s]

  upl1 <- apply(z21s, c(1, 2),
                function(z2) {
                  -sum((y - X%*%z2)^2)/(2*sigma^2) -
                    sum(abs(sqrt(gamma(3/q)/gamma(1/q))*z2/tau)^q)
                })
  upl2 <- apply(z22s, c(1, 2),
                function(z2) {
                  -sum((y - X%*%z2)^2)/(2*sigma^2) -
                    sum(abs(sqrt(gamma(3/q)/gamma(1/q))*z2/tau)^q)
                })
  upl3 <- apply(z23s, c(1, 2),
                function(z2) {
                  -sum((y - X%*%z2)^2)/(2*sigma^2) -
                    sum(abs(sqrt(gamma(3/q)/gamma(1/q))*z2/tau)^q)
                })

  dens1 <- lapply(1:chains, function(chain) {
    density(upl1[, chain])
  })
  dens2 <- lapply(1:chains, function(chain) {
    density(upl2[, chain])
  })
  dens3 <- lapply(1:chains, function(chain) {
    density(upl3[, chain])
  })

  xlim <- range(unlist(lapply(dens1, function(x) {x$x})),
                unlist(lapply(dens2, function(x) {x$x})),
                unlist(lapply(dens3, function(x) {x$x})))
  ylim <- range(unlist(lapply(dens1, function(x) {x$y})),
                unlist(lapply(dens2, function(x) {x$y})),
                unlist(lapply(dens3, function(x) {x$y})))

  div1 <- unlist(lapply(fitf1@sim[[1]], function(fit) {sum(attr(fit, "sampler_params")$divergent__[-(1:burn)])}))
  div2 <- unlist(lapply(fitf2@sim[[1]], function(fit) {sum(attr(fit, "sampler_params")$divergent__[-(1:burn)])}))
  div3 <- unlist(lapply(fitf3@sim[[1]], function(fit) {sum(attr(fit, "sampler_params")$divergent__[-(1:burn)])}))

  tim1 <- unlist(lapply(fitf1@sim[[1]], function(fit) {sum(attr(fit, "elapsed_time"))}))
  tim2 <- unlist(lapply(fitf2@sim[[1]], function(fit) {sum(attr(fit, "elapsed_time"))}))
  tim3 <- unlist(lapply(fitf3@sim[[1]], function(fit) {sum(attr(fit, "elapsed_time"))}))

  ess <- cbind(ess1, ess2, ess3)
  div <- cbind(div1, div2, div3)
  tim <- cbind(tim1, tim2, tim3)

  save(q, div, tim, ess, upl1, upl2, upl3, file = paste(dir, "Out/Data/simple_", dataname, "_q", i, ".RData", sep = ""))

}
