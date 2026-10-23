f0 <- function(z, mean0 = 0, sd0 = 1) {
  dnorm(z, mean = mean0, sd = sd0)
}

f.mixture <- function(z, pi0, mean1 = 0, sd1 = 1, mean0 = 0, sd0 = 1) { # nolint: object_name_linter
  pi0 * f0(z, mean0, sd0) + (1 - pi0) * dnorm(z, mean = mean1, sd = sd1)
}

fdr <- function(z, pi0, mean1, sd1, mean0 = 0, sd0 = 1) {
  # local false discovery rate for a test statistic = z
  # pi0: number of null-features
  # mean1: mean of distribution under H1
  # sd1: standard deviation under H1
  # mean0: mean under H0
  # st0: standard deviation under H0
  pi0 * f0(z, mean0, sd0) / f.mixture(z, pi0, mean1, sd1, mean0, sd0) # nolint: object_name_linter
}

rbnb <- function(n, mu, phi, rho) {
  # function for generating a bivariate NB with corr=rho
  # based on Gaussian copula
  # mu: marginal mean
  # phi: overdispersion (the same for both dimensions)
  # rho: target Spearman rank correlation
  rho <- (6 / pi) * asin(rho / 2)
  cop <- copula::normalCopula(param = rho, dim = 2)
  U <- copula::rCopula(n, cop) # nolint: object_name_linter
  X1 <- qnbinom(U[, 1], size = 1 / phi, mu = mu[1]) # nolint: object_name_linter
  X2 <- qnbinom(U[, 2], size = 1 / phi, mu = mu[2]) # nolint: object_name_linter
  cbind(X1, X2)
}

pbnb <- function(q, size, mu, rho) {
  #' CDF of a bivariate NB via Gaussian copula with target Spearman rho_s
  #'
  #' @param q   =(q1,q2) integer quantile for X1 (P(X1 <= q1)) and X2
  #' @param phi common NB size parameter r (>0)
  #' @param mu  marginal means of X1 and X2 (>0)
  #' @param rho desired Spearman rank correlation (-1 < rho < 1)
  #' @return     approximate CDF P(X1 <= q1, X2 <= q2)
  q1 <- q[1] # nolint: object_name_linter
  q2 <- q[2] # nolint: object_name_linter
  rho0 <- 2 * sin(pi * rho / 6) # nolint: object_name_linter
  u1 <- pnbinom(q1, size = size, mu = mu[1]) # nolint: object_name_linter
  u2 <- pnbinom(q2, size = size, mu = mu[2]) # nolint: object_name_linter
  upper <- qnorm(c(u1, u2))
  Sigma <- matrix(c(1, rho0, rho0, 1), 2, 2) # nolint: object_name_linter
  as.numeric(mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf),
    upper = upper,
    sigma = Sigma
  ))
}

cov.uv <- function(n0, n1, auc.u, auc.v, rho, distr = "norm", pars = NULL) { # nolint: object_name_linter
  # calculates the covariance between two WMW test statistics T_u and T_v
  n <- n0 + n1
  p.u <- auc.u # nolint: object_name_linter
  p.v <- auc.v # nolint: object_name_linter
  delta.u <- qnorm(auc.u) * sqrt(2) # nolint: object_name_linter
  delta.v <- qnorm(auc.v) * sqrt(2) # nolint: object_name_linter
  mu.1 <- c(delta.u, delta.v) # nolint: object_name_linter
  mu.0 <- c(0, 0) # nolint: object_name_linter
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2) # nolint: object_name_linter

  if (distr == "t") {
    z0a <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
    z0b <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
    zdiff <- z0a - z0b # nolint: object_name_linter
    delta.u <- as.numeric(quantile(zdiff, auc.u)) # nolint: object_name_linter
    delta.v <- as.numeric(quantile(zdiff, auc.v)) # nolint: object_name_linter
    mu.1 <- c(delta.u, delta.v) # nolint: object_name_linter
  }

  if (distr == "nbinom") {
    z0a <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
    z0b <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
    zdiff <- z0a / (z0b + 0.5) # nolint: object_name_linter
    delta.u <- as.numeric(quantile(zdiff, auc.u, na.rm = TRUE)) # nolint: object_name_linter
    delta.v <- as.numeric(quantile(zdiff, auc.v, na.rm = TRUE)) # nolint: object_name_linter
    mu.0 <- c(pars[1], pars[1]) # nolint: object_name_linter
    mu.1 <- mu.0 * c(delta.u, delta.v) # nolint: object_name_linter
  }

  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2) # nolint: object_name_linter

  uv0 <- switch(distr, # nolint: object_name_linter
    "norm"   = MASS::mvrnorm(10000, mu = mu.0, Sigma = Sigma),
    "t"      = mvtnorm::rmvt(10000, sigma = Sigma, df = pars),
    "nbinom" = rbnb(10000, mu = mu.0, phi = pars[2], rho = rho)
  )
  tau.uv.0 <- switch(distr, # nolint: object_name_linter
    "norm" = 1 - (1 - p.u) - (1 - p.v) +
      mean(pnorm(uv0[, 1], mean = mu.1[1]) *
             pnorm(uv0[, 2], mean = mu.1[2])) - p.u * p.v,
    "t" = 1 - (1 - p.u) - (1 - p.v) +
      mean(pt(uv0[, 1] - mu.1[1], df = pars) *
             pt(uv0[, 2] - mu.1[2], df = pars)) - p.u * p.v,
    "nbinom" = 1 - (1 - p.u) - (1 - p.v) +
      mean(pnbinom(uv0[, 1], mu = mu.1[1], size = 1 / pars[2]) *
             pnbinom(uv0[, 2], mu = mu.1[2], size = 1 / pars[2])) -
      p.u * p.v
  )

  uv1 <- switch(distr, # nolint: object_name_linter
    "norm"   = MASS::mvrnorm(10000, mu = mu.1, Sigma = Sigma),
    "t"      = mvtnorm::rmvt(10000, sigma = Sigma, df = pars,
      delta = mu.1, type = "shifted"), # nolint: indentation_linter
    "nbinom" = rbnb(10000, mu = mu.1, phi = pars[2], rho = rho)
  )
  tau.uv.1 <- switch(distr, # nolint: object_name_linter
    "norm" = mean(
      pnorm(uv1[, 1], mean = mu.0[1]) *
        pnorm(uv1[, 2], mean = mu.0[2])
    ) - p.u * p.v,
    "t" = mean(
      pt(uv1[, 1] - mu.0[1], df = pars) *
        pt(uv1[, 2] - mu.0[2], df = pars)
    ) - p.u * p.v,
    "nbinom" = mean(
      pnbinom(uv1[, 1], mu = mu.0[1], size = 1 / pars[2]) *
        pnbinom(uv1[, 2], mu = mu.0[2], size = 1 / pars[2])
    ) - p.u * p.v
  )

  Fuv0.eval <- switch(distr, # nolint: object_name_linter
    "norm"   = apply(uv1, 1, function(uv) {
      mvtnorm::pmvnorm(upper = uv, mean = mu.0, sigma = Sigma)
    }),
    "t"      = apply(uv1, 1, function(uv) {
      mvtnorm::pmvt(upper = uv, sigma = Sigma, df = pars)
    }),
    "nbinom" = apply(uv1, 1, function(uv) {
      pbnb(uv, size = 1 / pars[2], mu = mu.0, rho = rho)
    })
  )
  tau.uv <- mean(Fuv0.eval) - p.u * p.v # nolint: object_name_linter
  cov.uv <- 12 * ((n1 - 1) * tau.uv.0 + (n0 - 1) * tau.uv.1 + # nolint: object_name_linter
                    tau.uv) / (n + 1)
  cov.uv # nolint: object_name_linter
}


lfdr.sims <- function(scenario) { # nolint: object_name_linter
  # Simulates test statistics under H0 and H1 and calculates local fdrs.
  # All settings are specified in "scenario".
  # if scenario$u is NULL:
  #   if scenario$corr is NULL: test statistics are independently distributed
  #   if scenario$corr is not NULL: compound symmetry dependence structure
  # if scenario$u is not NULL: specifies an "arbitrary" dependence structure

  p  <- scenario$p
  p1 <- scenario$p1
  auc   <- scenario$auc
  alpha <- scenario$alpha # nolint: object_usage_linter
  n0    <- scenario$n0
  n1    <- scenario$n1
  auc0  <- scenario$auc0 # nolint: object_name_linter
  auc1  <- scenario$auc1 # nolint: object_name_linter
  B     <- scenario$B # nolint: object_name_linter
  u     <- scenario$u
  u.sd  <- scenario$u.sd # nolint: object_name_linter
  u.ind.fun <- ifelse( # nolint: object_name_linter
    is.null(scenario$u.ind.fun), NA, get(scenario$u.ind.fun)) # nolint: indentation_linter
  u.start.seq <- scenario$u.start.seq # nolint: object_name_linter
  pars  <- scenario$pars

  rho <- ifelse(is.null(scenario$corr), 0, scenario$corr)

  pi0 <- (p - p1) / p
  n   <- n0 + n1
  var0_w <- (n + 1) / (12 * n1 * n0) # nolint: object_name_linter

  if (is.na(auc0)) {
    if (scenario$distr == "norm") {
      delta <- sqrt(2) * qnorm(scenario$auc)
      z0a <- rnorm(1000000, mean = 0) # nolint: object_name_linter
      z0b <- rnorm(1000000, mean = 0) # nolint: object_name_linter
      z1a <- rnorm(1000000, mean = delta) # nolint: object_name_linter
      z1b <- rnorm(1000000, mean = delta) # nolint: object_name_linter
    }
    if (scenario$distr == "t") {
      delta <- find.delta.T(auc, pars) # nolint: object_name_linter
      z0a <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
      z0b <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
      zdiff <- z0a - z0b # nolint: object_name_linter
      delta <- as.numeric(quantile(zdiff, scenario$auc))
      z1a <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) + delta # nolint: object_name_linter
      z1b <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) + delta # nolint: object_name_linter
    }
    if (scenario$distr == "nbinom") {
      delta <- find.delta.NB(auc, pars) # nolint: object_name_linter
      z0a <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
      z0b <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
      z1a <- rnbinom(1000000, mu = pars[1] * delta, # nolint: object_name_linter
                     size = 1 / pars[2])
      z1b <- rnbinom(1000000, mu = pars[1] * delta, # nolint: object_name_linter
                     size = 1 / pars[2])
    }
    auc0 <- mean(pmax(z0a, z0b) < z1a) # nolint: object_name_linter
    auc1 <- mean(pmax(z1a, z1b) < z0a) # nolint: object_name_linter
  }

  var_w <- (((n0 - 1) * auc0 + (n1 - 1) * auc1) - # nolint: object_name_linter
              (n0 - 1) * auc^2 - (n1 - 1) * (1 - auc)^2 +
              auc * (1 - auc)) / (n0 * n1)
  sd1 <- sqrt(var_w) # nolint: object_name_linter
  sd0 <- sqrt(var0_w) # nolint: object_name_linter
  mean1 <- (auc - 0.5) / sd0 # nolint: object_name_linter
  var_w_upper  <- auc * (1 - auc) * max(n0, n1) / (n0 * n1) # nolint: object_name_linter
  sd1.upper    <- sqrt(var_w_upper) # nolint: object_name_linter

  local.fdr.null        <- matrix(nrow = B, ncol = p - p1) # nolint: object_name_linter
  local.fdr.pred        <- matrix(nrow = B, ncol = p1) # nolint: object_name_linter
  local.fdr.pred.upper  <- matrix(nrow = B, ncol = p1) # nolint: object_name_linter
  test.stats            <- matrix(nrow = B, ncol = p) # nolint: object_name_linter
  test.stats.upper      <- matrix(nrow = B, ncol = p) # nolint: object_name_linter

  if (is.null(u.sd)) {
    if (rho != 0) {
      tau.11 <- cov.uv( # nolint: object_name_linter
        n0 = n0, n1 = n1, auc.u = auc, auc.v = auc, rho = rho) # nolint: indentation_linter
      tau.00 <- cov.uv( # nolint: object_name_linter
        n0 = n0, n1 = n1, auc.u = .5, auc.v = .5, rho = rho) # nolint: indentation_linter
      tau.01 <- cov.uv( # nolint: object_name_linter
        n0 = n0, n1 = n1, auc.u = 0.5, auc.v = auc, rho = rho) # nolint: indentation_linter
      gamma <- tau.01 / tau.00
    } else {
      tau.00 <- 0 # nolint: object_name_linter
      tau.11 <- 0 # nolint: object_name_linter
      tau.01 <- 0 # nolint: object_name_linter
      gamma  <- 0
    }
    for (i in 1:B) {
      M  <- rnorm(1, sd = sqrt(tau.00)) # nolint: object_name_linter
      E0 <- rnorm(p - p1, sd = sqrt(1 - tau.00)) # nolint: object_name_linter
      E1 <- rnorm(p1, sd = sqrt(sd1^2 / sd0^2 - tau.11)) # nolint: object_name_linter
      z  <- mean1 + gamma * M + E1
      z0 <- M + E0 # nolint: object_name_linter
      local.fdr.pred[i, ] <- fdr( # nolint: object_name_linter
        z, pi0, mean1 = mean1, sd1 = sd1 / sd0, mean0 = 0, sd0 = 1) # nolint: indentation_linter
      local.fdr.null[i, ] <- fdr( # nolint: object_name_linter
        z0, pi0, mean1 = mean1, sd1 = sd1 / sd0, mean0 = 0, sd0 = 1) # nolint: indentation_linter
      test.stats[i, ] <- c(z, z0) # nolint: object_name_linter
      E1 <- rnorm(p1, sd = sqrt(sd1.upper^2 / sd0^2 - tau.11)) # nolint: object_name_linter
      z.upper <- mean1 + gamma * M + E1 # nolint: object_name_linter
      local.fdr.pred.upper[i, ] <- fdr( # nolint: object_name_linter
        z.upper, pi0, mean1 = mean1, sd1 = sd1.upper / sd0, # nolint: indentation_linter
        mean0 = 0, sd0 = 1)
      test.stats.upper[i, ] <- c(z.upper, z0) # nolint: object_name_linter
    }
  }

  if (!is.null(u.sd)) {
    u <- qnorm(c(seq(u.start.seq, 1 - u.start.seq, length.out = p)),
               mean = 0, sd = u.sd)
    u <- tanh(u)
    u <- sort(u)
    u.ind.fun <- ifelse( # nolint: object_name_linter
      is.null(scenario$u.ind.fun), NA, get(scenario$u.ind.fun)) # nolint: indentation_linter
    ind.z1 <- u.ind.fun(scenario = scenario) # nolint: object_name_linter
    u <- c(u[ind.z1], u[-ind.z1])
    zeta00 <- 0.8802113851 + 0.0050156840 * auc + # nolint: object_name_linter
      0.0066406274 * n0 - 0.0001320253 * n0^2
    zeta01 <- -0.1749915129 + 3.9685194231 * auc + # nolint: object_name_linter
      0.0081506108 * n0 - 0.0001314472 * n0^2 -
      3.7727232665 * auc^2 - 0.0022553820 * auc * n0
    zeta11 <- 9.959209e-01 + 1.070431 * auc + # nolint: object_name_linter
      (9.447239e-03) * n0 - (9.745422e-05) * n0^2 -
      2.240246 * auc^2 - (6.200959e-03) * auc * n0
    gamma1 <- sqrt(zeta11) # nolint: object_name_linter
    gamma0 <- sqrt(zeta00) # nolint: object_name_linter
    gamma  <- c(rep(gamma1, p1), rep(gamma0, p - p1))
    u2     <- u^2 # nolint: object_name_linter
    wtilde <- gamma * u2 / (1 - u2) # nolint: object_name_linter
    vk     <- cumsum(wtilde) # nolint: object_name_linter
    vkt    <- cumsum(gamma * wtilde) # nolint: object_name_linter
    var.t  <- 1 - u2[-1] * gamma[-1]^2 * # nolint: object_name_linter
      (vkt[-p] - (vk[-p]^2) / (1 + vkt[-p]))
    scenario$u <- u

    for (i in 1:B) {
      t.sim <- rnorm(1, mean = 0, sd = 1) # nolint: object_name_linter
      assign("t.sim", value = t.sim, pos = 1)
      sapply(2:p, FUN = function(j) {
        tmp <- wtilde[1:(j - 1)] * t.sim / u[1:(j - 1)]
        t.new <- rnorm(1, # nolint: object_name_linter
          mean = u[j] * gamma[j] / (1 + vkt[j - 1]) *
            sum(tmp[!is.na(tmp)]),
          sd = sqrt(var.t[(j - 1)])
        )
        assign("t.sim", value = c(get("t.sim", pos = 1), t.new), pos = 1)
      })
      t.sim <- get("t.sim", pos = 1) # nolint: object_name_linter
      z  <- t.sim[1:p1] * sd1 / sd0 + mean1
      z0 <- t.sim[-(1:p1)] # nolint: object_name_linter
      local.fdr.pred[i, ] <- fdr( # nolint: object_name_linter
        z, pi0, mean1 = mean1, sd1 = sd1 / sd0, mean0 = 0, sd0 = 1) # nolint: indentation_linter
      local.fdr.null[i, ] <- fdr( # nolint: object_name_linter
        z0, pi0, mean1 = mean1, sd1 = 1, mean0 = 0, sd0 = 1) # nolint: indentation_linter
      test.stats[i, ] <- c(z, z0) # nolint: object_name_linter
      z.upper <- t.sim[1:p1] * sd1.upper / sd0 + mean1 # nolint: object_name_linter
      local.fdr.pred.upper[i, ] <- fdr( # nolint: object_name_linter
        z.upper, pi0, mean1 = mean1, sd1 = sd1 / sd0, mean0 = 0, sd0 = 1) # nolint: indentation_linter
      test.stats.upper[i, ] <- c(z, z0) # nolint: object_name_linter
    }
  }

  list(
    lfdr.pred        = local.fdr.pred, # nolint: object_name_linter
    lfdr.null        = local.fdr.null, # nolint: object_name_linter
    test.stats       = test.stats, # nolint: object_name_linter
    lfdr.pred.upper  = local.fdr.pred.upper, # nolint: object_name_linter
    test.stats.upper = test.stats.upper, # nolint: object_name_linter
    tau = list(tau.00 = tau.00, tau.11 = tau.11, tau.01 = tau.01) # nolint: object_name_linter
  )
}


sens <- function(lfdr, scenario) {
  # Computes sensitivity from simulation results of lfdr.sims
  # lfdr: vector resulting from lfdr.sim

  alpha    <- scenario$alpha
  p        <- scenario$p
  p1       <- scenario$p1
  fdr.type <- scenario$fdr.type # nolint: object_name_linter

  if (fdr.type == "FDR") { # nolint: object_name_linter
    fd <- cbind(lfdr$lfdr.pred, lfdr$lfdr.null)
    tpp <- sapply(seq_len(nrow(fd)), function(i) {
      positive <- rep(FALSE, p)
      Fdr <- (cumsum(fd[i, order(lfdr$test.stats[i, ], # nolint: object_name_linter,indentation_linter
                                  decreasing = TRUE)]) / (1:p)) # nolint: indentation_linter
      positive[order(lfdr$test.stats[i, ], # nolint: object_name_linter
                     decreasing = TRUE)[Fdr < 0.05]] <- TRUE
      mean(positive[1:p1])
    })
    fpp <- sapply(seq_len(nrow(fd)), function(i) {
      positive <- rep(FALSE, p)
      Fdr <- (cumsum(fd[i, order(lfdr$test.stats[i, ], # nolint: object_name_linter
                                  decreasing = TRUE)]) / (1:p)) # nolint: indentation_linter
      positive[order(lfdr$test.stats[i, ], # nolint: object_name_linter
                     decreasing = TRUE)[Fdr < 0.05]] <- TRUE
      FPP <- sum(positive[-(1:p1)]) / max(sum(positive), 1) # nolint: object_name_linter
      FPP # nolint: object_name_linter
    })
    tpp.upper <- sapply(seq_len(nrow(fd)), function(i) { # nolint: object_name_linter
      positive <- rep(FALSE, p)
      Fdr <- (cumsum(fd[i, order(lfdr$test.stats.upper[i, ], # nolint: object_name_linter
                                  decreasing = TRUE)]) / (1:p)) # nolint: indentation_linter
      positive[order(lfdr$test.stats.upper[i, ], # nolint: object_name_linter
                     decreasing = TRUE)[Fdr < 0.05]] <- TRUE
      mean(positive[1:p1])
    })
  }

  if (fdr.type == "fdr") { # nolint: object_name_linter
    tpp       <- apply(lfdr$lfdr.pred < alpha, 1, mean)
    tpp.upper <- apply(lfdr$lfdr.pred.upper < alpha, 1, mean) # nolint: object_name_linter
  }

  list(
    tpp       = mean(tpp),
    tpp.upper = mean(tpp.upper), # nolint: object_name_linter
    ptp1      = mean(tpp > 0),
    fpp       = mean(fpp)
  )
}


wmw.test <- function(y, a) { # nolint: object_name_linter
  # Calculates the standardised Wilcoxon-Mann-Whitney test statistic
  # under the strong null hypothesis of equal distributions
  # (quick implementation -- no ties correction yet)
  # y: vector with all outcomes
  # a: vector with 0/1 treatment indicators
  n1  <- sum(a)
  n0  <- length(y) - n1
  tst <- wilcox.test(y ~ a)$stat
  tst <- (tst - 0.5 * n0 * n1) / sqrt(n0 * n1 * (n0 + n1 + 1) / 12)
  tst
}

calc.auc01 <- function(scenario) { # nolint: object_name_linter
  # Computes AUC0 and AUC1 by simulation from a hypothesised distribution
  delta <- sqrt(2) * qnorm(scenario$auc)
  pars  <- scenario$pars
  if (scenario$distr == "norm") {
    z0a <- rnorm(1000000, mean = 0) # nolint: object_name_linter
    z0b <- rnorm(1000000, mean = 0) # nolint: object_name_linter
    z1a <- rnorm(1000000, mean = delta) # nolint: object_name_linter
    z1b <- rnorm(1000000, mean = delta) # nolint: object_name_linter
  }
  if (scenario$distr == "t") {
    z0a <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
    z0b <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
    zdiff <- z0a - z0b # nolint: object_name_linter
    delta <- as.numeric(quantile(zdiff, scenario$auc))
    z1a <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) + delta # nolint: object_name_linter
    z1b <- rt(1000000, df = pars) / sqrt(pars / (pars - 2)) + delta # nolint: object_name_linter
  }
  if (scenario$distr == "nbinom") {
    z0a <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
    z0b <- rnbinom(1000000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
    zdiff <- z0a / (z0b + 0.5) # nolint: object_name_linter
    delta <- as.numeric(quantile(zdiff, scenario$auc, na.rm = TRUE))
    auc.tmp <- 0 # nolint: object_name_linter
    err <- (auc.tmp - scenario$auc) # nolint: object_name_linter
    sign.swich <- FALSE # nolint: object_name_linter
    while ((abs(err) > 0.02) && (!sign.swich)) {
      sign.err <- sign(err) # nolint: object_name_linter
      delta <- delta - sign.err
      z1a <- rnbinom(1000000, mu = pars[1] * delta, # nolint: object_name_linter
                     size = 1 / pars[2])
      auc.tmp <- mean(z1a > z0a) # nolint: object_name_linter
      err <- (auc.tmp - scenario$auc) # nolint: object_name_linter
      sign.swich <- (sign(err) != sign.err) # nolint: object_name_linter
    }
    z1a <- rnbinom(1000000, mu = pars[1] * delta, size = 1 / pars[2]) # nolint: object_name_linter
    z1b <- rnbinom(1000000, mu = pars[1] * delta, size = 1 / pars[2]) # nolint: object_name_linter
  }
  auc0 <- mean(pmax(z0a, z0b) < z1a) # nolint: object_name_linter
  auc1 <- mean(pmax(z1a, z1b) < z0a) # nolint: object_name_linter
  list(auc0 = auc0, auc1 = auc1) # nolint: object_name_linter
}


study.size <- function(scenario, sens.target = 0.8) { # nolint: object_name_linter
  # Computes the number of features that gives the desired sensitivity.
  # scenario$p will be ignored.
  # sens.target: the required sensitivity

  if (is.na(scenario$auc0)) {
    tmp <- calc.auc01(scenario)
    scenario$auc0 <- tmp$auc0
    scenario$auc1 <- tmp$auc1
  }

  scenario$n1 <- scenario$n0
  if (!is.null(scenario$n1.factor)) {
    scenario$n1 <- round(scenario$n0 * scenario$n1.factor)
  }
  scenario$p <- scenario$p1
  fdr.type   <- scenario$fdr.type # nolint: object_name_linter

  max.p <- 500000 # nolint: object_name_linter

  p.all <- c( # nolint: object_name_linter
    scenario$p1, seq(10, 100, 10), 125, seq(150, 500, 50), 750,
    seq(1000, 5000, 500),
    seq(6000, 10000, 1000), 15000, seq(20000, 50000, 10000),
    75000,
    seq(100000, 500000, 50000), max.p + 1
  )
  p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter

  scenario$p <- scenario$p1 + 1
  lf <- lfdr.sims(scenario = scenario)
  s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
  if (s1 > sens.target) {
    scenario$p <- 100
    lf <- lfdr.sims(scenario = scenario)
    s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
    scenario$p <- scenario$p1 + 1
    p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter
  }
  if (s1 > sens.target) {
    scenario$p <- 1000
    lf <- lfdr.sims(scenario = scenario)
    s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
    scenario$p <- 100
    p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter
  }
  if (s1 > sens.target) {
    scenario$p <- 10000
    lf <- lfdr.sims(scenario = scenario)
    s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
    scenario$p <- 1000
    p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter
  }
  if (s1 > sens.target) {
    scenario$p <- 100000
    lf <- lfdr.sims(scenario = scenario)
    s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
    scenario$p <- 10000
    p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter
  }
  if (s1 > sens.target) {
    scenario$p <- 200000
    lf <- lfdr.sims(scenario = scenario)
    s1 <- sens(lfdr = lf, scenario = scenario)$tpp # nolint: object_name_linter
    scenario$p <- 100000
    p.all <- p.all[p.all >= scenario$p] # nolint: object_name_linter
  }
  if (s1 > sens.target) {
    p.all <- p.all[p.all >= 200000] # nolint: object_name_linter
  }

  start.time <- Sys.time() # nolint: object_name_linter
  too.long   <- FALSE # nolint: object_name_linter
  s  <- 1
  i  <- 1
  while (s > sens.target) {
    p.previous <- scenario$p # nolint: object_name_linter
    scenario$p <- p.all[i]
    i <- i + 1
    lf <- lfdr.sims(scenario = scenario)
    s  <- sens(lfdr = lf, scenario = scenario)$tpp
    if (scenario$p > max.p) {
      too.long <- TRUE # nolint: object_name_linter
      break
    }
    now.time <- Sys.time() # nolint: object_name_linter
    if (difftime(now.time, start.time, units = "mins") > 10) {
      too.long <- TRUE # nolint: object_name_linter
      break
    }
  }
  scenario$p <- ifelse(p.previous == 1, 2, p.previous)

  s <- 1
  if (too.long) s <- 0
  start.time <- Sys.time() # nolint: object_name_linter
  while (s > sens.target) {
    p.previous  <- scenario$p # nolint: object_name_linter
    s.previous  <- s # nolint: object_name_linter
    scenario$p  <- round(scenario$p * 1.05) + 1
    lf <- lfdr.sims(scenario = scenario)
    s  <- sens(lfdr = lf, scenario = scenario)$tpp
    if (scenario$p > max.p) {
      too.long <- TRUE # nolint: object_name_linter
      break
    }
    now.time <- Sys.time() # nolint: object_name_linter
    if (difftime(now.time, start.time, units = "mins") > 10) {
      too.long <- TRUE # nolint: object_name_linter
      break
    }
  }
  scenario$p <- ifelse(p.previous == 1, 2, p.previous)

  di        <- floor(log10(p.previous))
  step.size <- max(10^(di - 1), 1) # nolint: object_name_linter
  s <- 1
  if (too.long) s <- 0
  start.time <- Sys.time() # nolint: object_name_linter
  while (s > sens.target) {
    p.previous <- scenario$p # nolint: object_name_linter
    s.previous <- s # nolint: object_name_linter
    scenario$p <- scenario$p + step.size
    lf <- lfdr.sims(scenario = scenario)
    s  <- sens(lfdr = lf, scenario = scenario)$tpp
    if (scenario$p > max.p) {
      p.previous <- NULL # nolint: object_name_linter
      break
    }
    now.time <- Sys.time() # nolint: object_name_linter
    if (difftime(now.time, start.time, units = "mins") > 10) {
      too.long <- TRUE # nolint: object_name_linter
      break
    }
  }
  scenario$p <- p.previous

  if (!too.long) {
    lf  <- lfdr.sims(scenario = scenario)
    tmp <- sens(lfdr = lf, scenario = scenario)
    sens.real  <- tmp$tpp # nolint: object_name_linter
    PTP1.real  <- tmp$ptp1 # nolint: object_name_linter
    FDR.real   <- tmp$fpp # nolint: object_name_linter
    res <- list(
      p          = scenario$p,
      sens.real  = sens.real, # nolint: object_name_linter
      PTP1.real  = PTP1.real, # nolint: object_name_linter
      FDR.real   = FDR.real # nolint: object_name_linter
    )
  }
  if (too.long) {
    res <- list(p = NA, sens.real = NA, PTP1.real = NA, FDR.real = NA) # nolint: object_name_linter
  }
  res
}


study.size.seq <- function(scenario, scenario.seq, # nolint: object_name_linter
                            sens.target = 0.8, path = NULL, graph = TRUE) { # nolint: indentation_linter
  # Sample size calculation for a sequence of scenarios.
  # scenario: list with scenario settings (except p1, n0 and auc)
  # scenario.seq: list with sequence for p1, n0 and auc
  # path: if a valid path is given, results will be saved
  # graph: if TRUE a graph will be produced
  p1.seq  <- scenario.seq$p1 # nolint: object_name_linter
  n0.seq  <- scenario.seq$n0 # nolint: object_name_linter
  auc.seq <- scenario.seq$auc.seq # nolint: object_name_linter

  fname.tmp <- tempfile( # nolint: object_name_linter
    pattern = "SSFSResults", # nolint: indentation_linter
    tmpdir = paste(path, "/SampleSizeResults", sep = ""))
  fname0 <- paste(fname.tmp, "_progress.txt", sep = "") # nolint: object_name_linter
  fname  <- paste(fname.tmp, ".RData", sep = "") # nolint: object_name_linter

  pars <- expand.grid(p1.seq, auc.seq, n0.seq)
  names(pars) <- c("p1", "auc", "n0")

  n.sc   <- nrow(pars) # nolint: object_name_linter
  n.iter <- floor(n.sc / ncores) # nolint: object_name_linter
  n.extra <- n.sc - n.iter * ncores # nolint: object_name_linter
  strt   <- 1
  np.results <- matrix(nrow = 1, ncol = 7) # nolint: object_name_linter
  if (n.iter != 0) {
    for (j in 1:n.iter) {
      results <- parallel::mclapply(strt:(strt + ncores - 1), # nolint: object_usage_linter
        function(i) {
          scenario$p1  <- pars$p1[i]
          scenario$n0  <- pars$n0[i]
          scenario$auc <- pars$auc[i]
          res <- study.size(scenario, sens.target = sens.target) # nolint: object_name_linter
          c(pars[i, ], res$p, res$sens.real, res$PTP1.real, res$FDR.real)
        },
        mc.cores = ncores
      )
      tmp <- do.call(rbind, results)
      np.results <- rbind(np.results, tmp) # nolint: object_name_linter
      strt <- strt + ifelse(is.null(nrow(tmp)), 0, nrow(tmp))
      if (!is.null(path)) {
        res <- write.results( # nolint: object_name_linter
          np.results = np.results, scenario = scenario,
          scenario.seq = scenario.seq, sens.target = sens.target,
          fname = fname, path = path
        )
      }
    }
  }
  strt <- n.iter * ncores + 1 # nolint: object_usage_linter
  if (n.extra > 0) {
    for (i in strt:(strt + n.extra - 1)) {
      scenario$p1  <- pars$p1[i]
      scenario$n0  <- pars$n0[i]
      scenario$auc <- pars$auc[i]
      tmp <- as.numeric(c(pars[i, ],
                          study.size(scenario, sens.target = sens.target)))
      np.results <- rbind(np.results, tmp) # nolint: object_name_linter
      if (!is.null(path)) {
        res <- write.results( # nolint: object_name_linter
          np.results = np.results, scenario = scenario,
          scenario.seq = scenario.seq, sens.target = sens.target,
          fname = fname, path = path
        )
      }
    }
  }

  if (!is.null(path)) {
    res <- write.results( # nolint: object_name_linter
      np.results = np.results, scenario = scenario,
      scenario.seq = scenario.seq, sens.target = sens.target,
      fname = fname, path = path
    )
  }

  if (graph) {
    ggplot2::ggplot(
      res,
      ggplot2::aes(
        x = n0, y = p, # nolint: object_usage_linter
        group = interaction(auc, p1), # nolint: object_usage_linter
        color = factor(auc),
        linetype = factor(p1)
      )
    ) +
      ggplot2::scale_y_continuous(trans = "log10") +
      ggplot2::geom_line()
  }

  invisible(res)
}


write.results <- function(np.results, scenario, scenario.seq, # nolint: object_name_linter
                           sens.target, fname, path) { # nolint: indentation_linter
  np.results <- as.data.frame(np.results) # nolint: object_name_linter
  np.results <- np.results[-1, ] # nolint: object_name_linter
  names(np.results) <- c("p1", "auc", "n0", "p", "sens.real", # nolint: object_name_linter
                          "PTP1.real", "FDR.real") # nolint: indentation_linter
  np.results$p1        <- unlist(np.results$p1) # nolint: object_name_linter
  np.results$auc       <- unlist(np.results$auc) # nolint: object_name_linter
  np.results$n0        <- unlist(np.results$n0) # nolint: object_name_linter
  np.results$p         <- unlist(np.results$p) # nolint: object_name_linter
  np.results$sens.real <- unlist(np.results$sens.real) # nolint: object_name_linter
  np.results$PTP1.real <- unlist(np.results$PTP1.real) # nolint: object_name_linter
  np.results$FDR.real  <- unlist(np.results$FDR.real) # nolint: object_name_linter
  np.results$n1        <- unlist(np.results$n0) # nolint: object_name_linter
  if (!is.null(scenario$n1.factor)) {
    scenario$n1.factor   <- unlist(scenario$n1.factor)
    np.results$n1        <- round(np.results$n0 * scenario$n1.factor) # nolint: object_name_linter
  }

  SSFSResults <- list( # nolint: object_name_linter
    scenario     = scenario,
    scenario.seq = scenario.seq, # nolint: object_name_linter
    sens.target  = sens.target, # nolint: object_name_linter
    results      = np.results
  )

  save(SSFSResults, file = fname) # nolint: object_name_linter
  invisible(np.results)
}


plot.study.size <- function(np, auc = NULL, p1 = NULL, # nolint: object_name_linter
                             n.range = NULL, p.range = NULL, # nolint: indentation_linter
                             plot.title = NULL,
                             log.scale = TRUE,
                             sens1 = FALSE) {
  # Makes graphs of p versus n from study size results matrix.
  # np: matrix with columns: p1, auc, n0, p and sens.real
  # auc: AUCs to include. If NULL, all AUCs in np are plotted.
  # p1: p1s to include. If NULL, all p1s in np are plotted.
  # n.range: range of n on x-axis. If NULL, full range in np.
  # p.range: range of p on x-axis. If NULL, full range in np.
  # plot.title: title printed on top of the plot.
  # log.scale: if TRUE, y-axis is on the log-scale.
  # sens1: if TRUE, results with real sensitivity of 100% are not plotted.

  if (!is.null(auc)) np <- np[np$auc %in% auc, ]
  if (!is.null(p1))  np <- np[np$p1 %in% p1, ]
  if (sens1)         np <- np[np$sens.real != 1, ]

  p <- ggplot2::ggplot(
    np,
    ggplot2::aes(
      x = n0, y = p, # nolint: object_usage_linter
      group    = interaction(auc, p1), # nolint: object_usage_linter
      color    = factor(auc),
      linetype = factor(p1)
    )
  ) +
    ggplot2::xlab("n (sample size per group)") +
    ggplot2::ylab("p (number of features)") +
    ggplot2::geom_line()

  if (log.scale) p <- p + ggplot2::scale_y_continuous(trans = "log10")
  if (!is.null(plot.title)) {
    p <- p + ggplot2::ggtitle(plot.title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5, vjust = 0.5, face = "bold"
        )
      )
  }
  if (!is.null(n.range)) p <- p + ggplot2::xlim(n.range)
  if (!is.null(p.range)) p <- p + ggplot2::ylim(p.range)
  p
}


find.study.size.results <- function(scenario, sens.target = 0.8, # nolint: object_name_linter
                                     path, ...) { # nolint: indentation_linter
  # Files in "path" will be searched for results for the given scenario.
  scen1 <- list(
    alpha       = scenario$alpha,
    n1.factor   = scenario$n1.factor, # nolint: object_name_linter
    u.sd        = scenario$u.sd, # nolint: object_name_linter
    u.start.seq = scenario$u.start.seq, # nolint: object_name_linter
    u.ind.fun   = scenario$u.ind.fun, # nolint: object_name_linter
    distr       = scenario$distr,
    pars        = scenario$pars,
    description = scenario$description,
    sens.target = sens.target, # nolint: object_name_linter
    corr        = ifelse(is.null(scenario$corr), 0, scenario$corr)
  )

  all.files <- dir(paste(path, "/SampleSizeResults", sep = "")) # nolint: object_name_linter
  agree <- FALSE
  i <- 1
  while ((!agree) && (i <= length(all.files))) {
    fl <- all.files[i]
    load(paste(path, "/SampleSizeResults/", fl, sep = ""))
    scen2 <- list(
      alpha       = SSFSResults$scenario$alpha, # nolint: object_name_linter
      n1.factor   = SSFSResults$scenario$n1.factor, # nolint: object_name_linter
      u.sd        = SSFSResults$scenario$u.sd, # nolint: object_name_linter
      u.start.seq = SSFSResults$scenario$u.start.seq, # nolint: object_name_linter
      u.ind.fun   = SSFSResults$scenario$u.ind.fun, # nolint: object_name_linter
      distr       = SSFSResults$scenario$distr, # nolint: object_name_linter
      pars        = SSFSResults$scenario$pars, # nolint: object_name_linter
      description = SSFSResults$scenario$description, # nolint: object_name_linter
      sens.target = SSFSResults$sens.target, # nolint: object_name_linter
      corr        = ifelse( # nolint: object_name_linter
        is.null(SSFSResults$scenario$corr), 0, SSFSResults$scenario$corr) # nolint: object_name_linter
    )
    if (length(setdiff(scen1, scen2)) == 0) {
      results <- SSFSResults$results # nolint: object_name_linter
      agree <- TRUE
    } else {
      remove("SSFSResults")
    }
    i <- i + 1
  }

  if (agree) {
    results <- SSFSResults$results # nolint: object_name_linter
    plot.study.size(results, ...) # nolint: object_name_linter
    SSFSResults # nolint: object_name_linter
  } else {
    warning("Results not found. \n")
  }
}


list.study.size.scenarios <- function(path) { # nolint: object_name_linter
  # Files in "path" will be listed.
  all.files <- dir(paste(path, "/SampleSizeResults", sep = "")) # nolint: object_name_linter
  results <- data.frame(
    filename    = all.files,
    alpha       = NA,
    n1.factor   = NA, # nolint: object_name_linter
    distr       = "",
    descr       = "",
    sens.target = NA, # nolint: object_name_linter
    corr        = NA,
    u.sd        = NA, # nolint: object_name_linter
    u.start.seq = NA, # nolint: object_name_linter
    u.ind.fun   = "" # nolint: object_name_linter
  )

  cnt <- 1
  for (fl in all.files) {
    load(paste(path, "/SampleSizeResults/", fl, sep = ""))
    results$alpha[cnt]       <- SSFSResults$scenario$alpha # nolint: object_name_linter
    results$n1.factor[cnt]   <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$n1.factor), 1, # nolint: object_name_linter
      SSFSResults$scenario$n1.factor) # nolint: object_name_linter
    results$u.sd[cnt]        <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$u.sd), NA, # nolint: object_name_linter
      SSFSResults$scenario$u.sd) # nolint: object_name_linter
    results$u.start.seq[cnt] <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$u.start.seq), NA, # nolint: object_name_linter
      SSFSResults$scenario$u.start.seq) # nolint: object_name_linter
    results$u.ind.fun[cnt]   <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$u.ind.fun), NA, # nolint: object_name_linter
      SSFSResults$scenario$u.ind.fun) # nolint: object_name_linter
    results$distr[cnt]       <- SSFSResults$scenario$distr # nolint: object_name_linter
    results$descr[cnt]       <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$description), NA, # nolint: object_name_linter
      SSFSResults$scenario$description) # nolint: object_name_linter
    results$sens.target[cnt] <- SSFSResults$sens.target # nolint: object_name_linter
    results$corr[cnt]        <- ifelse( # nolint: object_name_linter
      is.null(SSFSResults$scenario$corr), 0, # nolint: object_name_linter
      SSFSResults$scenario$corr) # nolint: object_name_linter
    cnt <- cnt + 1
  }
  results
}


### functions for general dependence

gen.dep.v1 <- function(cor.target = 0.3, scenario) { # nolint: object_name_linter
  u  <- scenario$u
  p1 <- scenario$p1
  u1 <- which.min((abs(u) - sqrt(cor.target))^2) # nolint: object_name_linter
  z.ind <- u1:(u1 + p1 - 1) # nolint: object_name_linter
  if (z.ind[p1] > scenario$p) {
    z.ind <- z.ind - (z.ind[p1] - scenario$p) # nolint: object_name_linter
  }
  z.ind # nolint: object_name_linter
}

gen.dep.v2 <- function(cor.target = 0.3, scenario) { # nolint: object_name_linter
  z.ind <- sample(1:scenario$p, scenario$p1, replace = FALSE) # nolint: object_name_linter
  z.ind # nolint: object_name_linter
}


move.doubles <- function(path) { # nolint: object_name_linter
  # Reads .RData files, identifies files with identical scenarios,
  # keeps the most recent one, and moves older duplicates.

  input_directory  <- paste(path, "/SampleSizeResults", sep = "") # nolint: object_name_linter
  output_directory <- paste(path, "/SampleSizeResultsDoubles", sep = "") # nolint: object_name_linter

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message(paste("Created output directory:", output_directory))
  }

  rdata_files <- list.files(input_directory, pattern = "\\.RData$", # nolint: object_name_linter
                             full.names = TRUE) # nolint: indentation_linter
  if (length(rdata_files) == 0) {
    stop("No .RData files found in the specified input directory.")
  }
  message(paste("Found", length(rdata_files), ".RData files."))

  scenarios_data <- list() # nolint: object_name_linter
  for (file_path in rdata_files) { # nolint: object_name_linter
    temp_env <- new.env() # nolint: object_name_linter
    load(file_path, envir = temp_env)
    if ("SSFSResults" %in% ls(envir = temp_env)) {
      scenario_content_tmp <- temp_env$SSFSResults$scenario # nolint: object_name_linter
      scenario_content <- list() # nolint: object_name_linter
      scenario_content$alpha       <- scenario_content_tmp$alpha # nolint: object_name_linter
      scenario_content$fdr.type    <- scenario_content_tmp$fdr.type # nolint: object_name_linter
      scenario_content$n0          <- scenario_content_tmp$n0 # nolint: object_name_linter
      scenario_content$distr       <- scenario_content_tmp$distr # nolint: object_name_linter
      scenario_content$pars        <- scenario_content_tmp$pars # nolint: object_name_linter
      scenario_content$description <- scenario_content_tmp$description # nolint: object_name_linter
      scenario_content$n1.factor   <- scenario_content_tmp$n1.factor # nolint: object_name_linter
      if (is.null(scenario_content$n1.factor)) # nolint: object_name_linter
        scenario_content$n1.factor <- 1 # nolint: object_name_linter
      if (is.na(scenario_content$n1.factor)) # nolint: object_name_linter
        scenario_content$n1.factor <- 1 # nolint: object_name_linter
      scenario_content$corr        <- scenario_content_tmp$corr # nolint: object_name_linter
      if (is.null(scenario_content$corr)) # nolint: object_name_linter
        scenario_content$corr <- 0 # nolint: object_name_linter
      if (is.na(scenario_content$corr)) # nolint: object_name_linter
        scenario_content$corr <- 0 # nolint: object_name_linter
      scenario_content$sens.target <- temp_env$SSFSResults$sens.target # nolint: object_name_linter
      scenarios_data[[file_path]] <- list( # nolint: object_name_linter
        scenario = scenario_content, # nolint: object_name_linter
        mtime    = file.info(file_path)$mtime
      )
      message(paste("Successfully loaded scenario from:",
                    basename(file_path)))
    } else {
      warning(paste("File", basename(file_path),
                    "does not contain 'scenario'. Skipping."))
    }
  }

  if (length(scenarios_data) == 0) {
    stop("No files contained a 'scenario' object. Exiting.")
  }

  scenario_info_df <- data.frame( # nolint: object_name_linter
    file_path   = character(), # nolint: object_name_linter
    mtime       = as.POSIXct(character()),
    scenario_id = character(), # nolint: object_name_linter
    stringsAsFactors = FALSE
  )

  for (file_path in names(scenarios_data)) { # nolint: object_name_linter
    scenario_content <- scenarios_data[[file_path]]$scenario # nolint: object_name_linter
    mtime            <- scenarios_data[[file_path]]$mtime
    scenario_json <- toJSON( # nolint: object_name_linter
      scenario_content, auto_unbox = TRUE, # nolint: indentation_linter
      prettify = FALSE, sort_by_name = TRUE)
    scenario_info_df <- rbind(scenario_info_df, data.frame( # nolint: object_name_linter
      file_path        = file_path, # nolint: object_name_linter
      mtime            = mtime,
      scenario_id      = scenario_json, # nolint: object_name_linter
      stringsAsFactors = FALSE
    ))
  }

  duplicate_scenario_ids <- scenario_info_df$scenario_id[ # nolint: object_name_linter
    duplicated(scenario_info_df$scenario_id)] # nolint: indentation_linter

  if (length(duplicate_scenario_ids) == 0) {
    message("No duplicate scenarios found.")
  } else {
    message(paste("Found", length(duplicate_scenario_ids),
                  "unique scenario types with duplicates."))
  }

  files_to_move <- c() # nolint: object_name_linter
  for (dup_id in unique(duplicate_scenario_ids)) { # nolint: object_name_linter
    matching_files_df <- scenario_info_df[ # nolint: object_name_linter
      scenario_info_df$scenario_id == dup_id, ] # nolint: indentation_linter
    matching_files_df <- matching_files_df[ # nolint: object_name_linter
      order(matching_files_df$mtime, decreasing = TRUE), ] # nolint: indentation_linter
    file_to_keep    <- matching_files_df$file_path[1] # nolint: object_name_linter
    message(paste("Keeping:", basename(file_to_keep)))
    older_duplicates <- matching_files_df$file_path[-1] # nolint: object_name_linter
    if (length(older_duplicates) > 0) {
      message(paste("Identified", length(older_duplicates),
                    "older duplicate(s)."))
      files_to_move <- c(files_to_move, older_duplicates) # nolint: object_name_linter
    }
  }

  if (length(files_to_move) == 0) {
    message("No files need to be moved.")
  } else {
    message(paste("Attempting to move", length(files_to_move), "files."))
    for (file_to_move in files_to_move) { # nolint: object_name_linter
      destination_path <- file.path(output_directory, # nolint: object_name_linter
                                    basename(file_to_move))
      if (file.exists(destination_path)) {
        warning(paste("File", basename(file_to_move),
                      "already exists in destination. Skipping."))
      } else {
        if (file.rename(file_to_move, destination_path)) {
          message(paste("Moved:", basename(file_to_move)))
        } else {
          warning(paste("Failed to move:", basename(file_to_move)))
        }
      }
    }
    message("File management complete.")
  }
}
