# basic simulation code.

sim.data <- function(n0 = 30, n1 = 30, p = 1000, p1 = 10, auc = 0.5, # nolint: object_name_linter
                     corr = NULL, u = NULL, cov.matrix = NULL, # nolint: object_name_linter
                     cov.matrix.sqrt = NULL, # nolint: object_name_linter
                     distr = "norm", pars = NA, delta = NA) {
  # Simulates a single dataset where the first p1 features are predictive.
  #
  # If cov.matrix = corr = NULL: features are independently distributed.
  # If corr != NULL: compound symmetry dependence structure.
  # If cov.matrix != NULL: arbitrary dependence structure
  #   (cov.matrix.sqrt must also be given).
  # distr: distribution to simulate from.
  # pars: vector with parameters for distr (t-distr: df; nbinom: c(mu, phi)).
  # delta: determined by AUC; for nbinom must be given explicitly.

  p0 <- p - p1 # nolint: object_name_linter

  db <- matrix(nrow = n0 + n1, ncol = p + 1)
  db <- as.data.frame(db)
  db[, -1] <- rnorm((n0 + n1) * p)

  if (!is.null(cov.matrix)) {
    db[, -1] <- as.matrix(db[, -1]) %*% cov.matrix.sqrt
  }
  if (!is.null(corr)) {
    e <- rnorm(n0 + n1, sd = sqrt(corr))
    db[, -1] <- db[, -1] * sqrt(1 - corr) + e
  }

  if (distr == "norm") {
    delta <- qnorm(auc) * sqrt(2)
  }

  db[, 1] <- c(rep(0, n0), rep(1, n1))

  if (distr == "norm") {
    db[1:n0, 2:(p1 + 1)] <- db[1:n0, 2:(p1 + 1)] + delta
  }
  if (distr == "t") {
    db[, -1] <- rt(p * (n0 + n1), df = pars) / sqrt(pars / (pars - 2))
    db[1:n0, 2:(p1 + 1)] <- db[1:n0, 2:(p1 + 1)] + delta
  }
  if (distr == "nbinom") {
    db[, -1] <- rnbinom(p * (n0 + n1), mu = pars[1], size = 1 / pars[2])
    db[1:n0, 2:(p1 + 1)] <- rnbinom(
      p1 * n0, mu = pars[1] * delta, size = 1 / pars[2]) # nolint: indentation_linter
  }

  names(db) <- c("a", paste("F", 1:p, sep = ""))
  db
}


deriv.auc <- function(x, derivs.at.midpoints) { # nolint: object_name_linter
  row <- which.min(abs(x - derivs.at.midpoints$auc))
  derivs.at.midpoints$d[row] # nolint: object_name_linter
}

sim1 <- function(scenario) {
  # Performs a simulation experiment for one single scenario.
  #
  # scenario$N: number of Monte Carlo runs
  # scenario$B: number of simulations for the lfdr.sim function
  #
  # Returns sensitivity, FDR, and bias/MSE results for
  # nonparametric and Tweedie estimators of the AUC.

  go.back <- FALSE # nolint: object_name_linter

  N        <- scenario$N # nolint: object_name_linter
  B        <- scenario$B # nolint: object_name_linter
  p        <- scenario$p
  p1       <- scenario$p1
  auc      <- scenario$auc
  alpha    <- scenario$alpha
  n0       <- scenario$n0
  n1       <- scenario$n1
  corr     <- scenario$corr
  fdr.type <- scenario$fdr.type # nolint: object_name_linter
  pars     <- scenario$pars
  delta    <- scenario$delta
  distr    <- scenario$distr

  test.stats.all    <- matrix(nrow = N, ncol = p) # nolint: object_name_linter
  lfdr.obs.all      <- matrix(nrow = N, ncol = p) # nolint: object_name_linter
  auc.all           <- matrix(nrow = N, ncol = p) # nolint: object_name_linter
  auc.Tweedie.all   <- matrix(nrow = N, ncol = p) # nolint: object_name_linter

  TPP <- FDP <- PTP1 <- numeric(N) # nolint: object_name_linter
  TPP.BH <- FDP.BH <- numeric(N) # nolint: object_name_linter
  for (i in 1:N) {
    db <- sim.data( # nolint: object_name_linter
      n0 = n0, n1 = n1, p = p, p1 = p1, auc = auc, corr = corr,
      u = scenario$u, cov.matrix = scenario$cov.matrix,
      cov.matrix.sqrt = scenario$cov.matrix.sqrt,
      distr = distr, pars = pars, delta = delta
    )

    test.stats <- apply(db[, -1], 2, wmw.test, a = db[, 1]) # nolint: object_name_linter
    test.stats.all[i, ] <- test.stats # nolint: object_name_linter

    p.adj.BH <- p.adjust(1 - pnorm(test.stats), method = "BH") # nolint: object_name_linter
    positive <- rep(FALSE, p)
    positive[p.adj.BH < alpha] <- TRUE
    TPP.BH[i] <- mean(positive[1:p1]) # nolint: object_name_linter
    FDP.BH[i] <- sum(positive[-(1:p1)]) / max(sum(positive), 1) # nolint: object_name_linter

    lfdr.obs <- try( # nolint: object_name_linter
      locfdr::locfdr(zz = test.stats, nulltype = 0, plot = 0, df = 8), # nolint: indentation_linter
      silent = TRUE)
    if (!inherits(lfdr.obs, "try-error")) {
      lfdr.obs.all[i, ] <- lfdr.obs$fdr # nolint: object_name_linter
      positive <- rep(FALSE, p)
      Fdr <- (cumsum(lfdr.obs$fdr[order(test.stats, # nolint: object_name_linter
                                         decreasing = TRUE)]) / (1:p)) # nolint: indentation_linter
      positive[order(test.stats,
                     decreasing = TRUE)[Fdr < alpha]] <- TRUE
      TPP[i]  <- mean(positive[1:p1]) # nolint: object_name_linter
      FDP[i]  <- sum(positive[-(1:p1)]) / max(sum(positive), 1) # nolint: object_name_linter
      PTP1[i] <- max(positive[1:p1]) # nolint: object_name_linter
    } else {
      go.back <- TRUE # nolint: object_name_linter
    }

    auc.tmp <- sapply(1:p, FUN = function(i, db) { # nolint: object_name_linter
      y  <- db[, i + 1] # nolint: object_name_linter
      a  <- db[, 1] # nolint: object_name_linter
      pm <- pim(y ~ a, link = "identity") # nolint: object_name_linter
      c(1 - pm@coef, pm@vcov)
    }, db = db)
    auc.all[i, ] <- auc.tmp[1, ] # nolint: object_name_linter

    auc.obs <- try( # nolint: object_name_linter
      locfdr::locfdr(zz = auc.tmp[1, ], nulltype = 0, plot = 0, df = 8), # nolint: indentation_linter
      silent = TRUE)
    if (!inherits(auc.obs, "try-error")) {
      mat <- auc.obs$mat
      dx  <- diff(mat[, 1])[1]
      derivs <- data.frame(auc = mat[-1, 1], d = diff(log(mat[, 5])) / dx)
      d <- sapply(auc.tmp[1, ], FUN = deriv.auc, # nolint: object_name_linter
                  derivs.at.midpoints = derivs)
      auc.Tweedie.all[i, ] <- auc.tmp[1, ] + auc.tmp[2, ] * d # nolint: object_name_linter
    } else {
      go.back <- TRUE # nolint: object_name_linter
    }

    if (go.back) {
      i <- i - 1 # nolint: object_name_linter
      go.back <- FALSE # nolint: object_name_linter
    }
  }

  AUC.estimation <- data.frame( # nolint: object_name_linter
    bias0          = mean(auc.all[, -(1:p1)] - 0.5),
    bias1          = mean(auc.all[, 1:p1] - auc),
    bias0.Tweedie  = mean(auc.Tweedie.all[, -(1:p1)] - 0.5), # nolint: object_name_linter
    bias1.Tweedie  = mean(auc.Tweedie.all[, 1:p1] - auc), # nolint: object_name_linter
    MSE0           = mean((auc.all[, -(1:p1)] - 0.5)^2),
    MSE1           = mean((auc.all[, 1:p1] - auc)^2),
    MSE0.Tweedie   = mean((auc.Tweedie.all[, -(1:p1)] - 0.5)^2), # nolint: object_name_linter
    MSE1.Tweedie   = mean((auc.Tweedie.all[, 1:p1] - auc)^2) # nolint: object_name_linter
  )

  sens.emp     <- mean(TPP) # nolint: object_name_linter
  FDR.emp      <- mean(FDP) # nolint: object_name_linter
  sens.emp.BH  <- mean(TPP.BH) # nolint: object_name_linter
  FDR.emp.BH   <- mean(FDP.BH) # nolint: object_name_linter
  PTP1.emp     <- mean(PTP1) # nolint: object_name_linter

  list(
    sens            = sens.emp, # nolint: object_name_linter
    FDR             = FDR.emp, # nolint: object_name_linter
    sens.BH         = sens.emp.BH, # nolint: object_name_linter
    FDR.BH          = FDR.emp.BH, # nolint: object_name_linter
    PTP1            = PTP1.emp, # nolint: object_name_linter
    test.stat.all   = test.stats.all, # nolint: object_name_linter
    lfdr.all        = lfdr.obs.all, # nolint: object_name_linter
    AUC.estimation  = AUC.estimation # nolint: object_name_linter
  )
}


sim.seq <- function(scenario, auc.seq = seq(0.5, 0.9, by = 0.1), # nolint: object_name_linter
                    p1.seq = seq(10, 110, by = 25), # nolint: object_name_linter
                    graph = TRUE, path = NULL) {
  # Performs a simulation study for a sequence of AUCs and p1s.
  # graph: should results be shown in a graph?
  # path: if a valid path is given, results will be saved.

  results <- as.data.frame(matrix(
    nrow = length(auc.seq) * length(p1.seq), ncol = 11
  ))
  names(results) <- c("auc", "p1", "sens.emp", "sens.new", "sens.new.upper",
                      "sens.emp.BH", "FDR", "FDR.BH", "MSE",
                      "MSE.Tweedie", "PTP1.emp")
  ncores        <- 6
  scenario$N    <- round(scenario$N / ncores)
  distr         <- scenario$distr
  u.start.seq   <- scenario$u.start.seq # nolint: object_name_linter
  u.sd          <- scenario$u.sd # nolint: object_name_linter
  p             <- scenario$p
  if (!is.null(scenario$u.sd)) {
    u <- qnorm(c(seq(u.start.seq, 1 - u.start.seq, length.out = p)),
               mean = 0, sd = u.sd)
    u <- tanh(u)
    u <- sort(u)
    u.ind.fun <- ifelse(is.null(scenario$u.ind.fun), NA, # nolint: object_name_linter
                        get(scenario$u.ind.fun))
    ind <- u.ind.fun(u, scenario = scenario)
    u   <- c(u[ind], u[-ind])
    U   <- u %*% t(u) # nolint: object_name_linter
    D   <- diag(1 - diag(U)) # nolint: object_name_linter
    scenario$cov.matrix      <- U + D
    scenario$cov.matrix.sqrt <- expm::sqrtm(scenario$cov.matrix)
  }
  cnt <- 1
  for (i.a in auc.seq) {
    if (distr == "nbinom") {
      scenario$delta <- find.delta.NB(auc = i.a, pars = scenario$pars) # nolint: object_name_linter
    }
    if (distr == "t") {
      scenario$delta <- find.delta.T(auc = i.a, pars = scenario$pars) # nolint: object_name_linter
    }
    for (i.p1 in p1.seq) {
      scenario$auc <- i.a
      scenario$p1  <- i.p1

      all.results <- parallel::mclapply(1:ncores, # nolint: object_name_linter
        function(i) {
          res   <- sim1(scenario)
          estim <- res$AUC.estimation
          c(
            res$sens, res$FDR,
            res$sens.BH, res$FDR.BH,
            (estim$MSE0 * (scenario$p - i.p1) +
               estim$MSE1 * i.p1) / scenario$p,
            (estim$MSE0.Tweedie * (scenario$p - i.p1) +
               estim$MSE1.Tweedie * i.p1) / scenario$p,
            res$PTP1
          )
        },
        mc.cores = ncores
      )

      lf       <- lfdr.sims(scenario = scenario) # nolint: object_name_linter
      sens.tmp <- sens(lfdr = lf, scenario = scenario) # nolint: object_name_linter
      sens.new       <- sens.tmp$tpp # nolint: object_name_linter
      sens.new.upper <- sens.tmp$tpp.upper # nolint: object_name_linter
      res <- do.call(rbind, all.results)
      results[cnt, ] <- c(
        i.a, i.p1,
        mean(res[, 1]), sens.new, sens.new.upper, mean(res[, 3]),
        mean(res[, 2]), mean(res[, 4]),
        mean(res[, 5]), mean(res[, 6]),
        mean(res[, 7])
      )
      cnt <- cnt + 1
      print(c(i.a, i.p1))
    }
  }

  if (graph) {
    plot(results$sens.emp, results$sens.new,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "empirical sensitivity",
      ylab = "computed sensitivity",
      cex.axis = 1.2, cex.lab = 1.2
    )
    abline(c(0, 1))
    plot(results$sens.emp, results$sens.new.upper,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "empirical sensitivity",
      ylab = "computed sensitivity",
      cex.axis = 1.2, cex.lab = 1.2,
      main = "upper bound of variance"
    )
    abline(c(0, 1))
    plot(results$sens.emp.BH, results$sens.new,
      xlim = c(0, 1), ylim = c(0, 1),
      xlab = "empirical sensitivity (BH)",
      ylab = "computed sensitivity",
      cex.axis = 1.2, cex.lab = 1.2,
      main = "BH adjusted"
    )
    abline(c(0, 1))
    mn <- min(results$MSE, results$MSE.Tweedie)
    mx <- max(results$MSE, results$MSE.Tweedie)
    plot(results$MSE, results$MSE.Tweedie,
      xlim = c(mn, mx), ylim = c(mn, mx),
      xlab = "MSE of the nonparametric AUC estimator",
      ylab = "MSE of the Tweedie AUC estimator",
      cex.axis = 1.2, cex.lab = 1.2
    )
    abline(c(0, 1))
  }

  if (!is.null(path)) {
    SimResults <- list(scenario = scenario, results = results) # nolint: object_name_linter
    save(SimResults, file = paste( # nolint: object_name_linter
      tempfile(pattern = "SimResults",
               tmpdir = paste(path, "/SimResults", sep = "")),
      ".RData", sep = ""
    ))
  }

  results
}


find.delta.NB <- function(auc, pars) { # nolint: object_name_linter
  # Computes the delta needed to reach auc (for the NB distribution).
  # auc: the desired AUC
  # pars: the parameters (mu, phi) for the NB

  x1     <- rnbinom(100000, mu = pars[1], size = 1 / pars[2]) # nolint: object_name_linter
  cnt    <- 1
  d      <- 1
  target <- FALSE
  while ((cnt < 1000) && (!target)) {
    d  <- d + 0.5
    x2 <- rnbinom(100000, mu = pars[1] * d, size = 1 / pars[2]) # nolint: object_name_linter
    auc2 <- mean(x1 < x2) + 0.5 * mean(x1 == x2) # nolint: object_name_linter
    if (auc2 >= auc) target <- TRUE
    cnt <- cnt + 1
  }
  if (cnt == 1000) stop("Could not reach the desired AUC \n")
  d
}


find.delta.T <- function(auc, pars) { # nolint: object_name_linter
  # Computes the delta needed to reach auc (for the T distribution).
  # auc: the desired AUC
  # pars: the parameter (df) for the T distribution

  x1     <- rt(100000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
  x2     <- rt(100000, df = pars) / sqrt(pars / (pars - 2)) # nolint: object_name_linter
  cnt    <- 1
  d      <- 0
  target <- FALSE
  while ((cnt < 5000) && (!target)) {
    d    <- d + 0.01
    auc2 <- mean(x1 < x2 + d) + 0.5 * mean(x1 == (x2 + d)) # nolint: object_name_linter
    if (auc2 >= auc) target <- TRUE
    cnt <- cnt + 1
  }
  if (cnt == 5000) stop("Could not reach the desired AUC \n")
  d
}


plot.sim.seq <- function(results) { # nolint: object_name_linter
  # Plots the results obtained with the sim.seq function.
  # results: output of the sim.seq function

  par(mfrow = c(2, 2))

  plot(results$sens.emp, results$sens.new,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "empirical sensitivity",
    ylab = "computed sensitivity",
    cex.axis = 1.2, cex.lab = 1.2,
    main = "use of AUC0 and AUC1"
  )
  abline(c(0, 1))

  plot(results$sens.emp.BH, results$sens.new,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "empirical sensitivity (BH)",
    ylab = "computed sensitivity",
    cex.axis = 1.2, cex.lab = 1.2,
    main = "BH adjusted"
  )
  abline(c(0, 1))

  plot(results$sens.emp, results$sens.new.upper,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "empirical sensitivity",
    ylab = "computed sensitivity",
    cex.axis = 1.2, cex.lab = 1.2,
    main = "use of upper bound"
  )
  abline(c(0, 1))

  mn <- min(c(results$MSE, results$MSE.Tweedie), na.rm = TRUE)
  mx <- max(c(results$MSE, results$MSE.Tweedie), na.rm = TRUE)
  plot(results$MSE, results$MSE.Tweedie,
    xlim = c(mn, mx), ylim = c(mn, mx),
    xlab = "MSE of the nonpar. AUC est.",
    ylab = "MSE of the Tweedie AUC est.",
    cex.axis = 1.2, cex.lab = 1.2
  )
  abline(c(0, 1))

  par(mfrow = c(1, 1))
}

find.sim.seq.results <- function(scenario, path) { # nolint: object_name_linter
  # Files in "path" will be searched for results for the given scenario.

  all.files <- dir(paste(path, "/SimResults", sep = "")) # nolint: object_name_linter
  agree <- FALSE
  i <- 1
  while ((!agree) && (i <= length(all.files))) {
    fl <- all.files[i]
    load(paste(path, "/SimResults/", fl, sep = ""))
    if (length(setdiff(scenario, SimResults$scenario)) == 0) { # nolint: object_name_linter
      results <- SimResults$results # nolint: object_name_linter
      agree <- TRUE
    } else {
      remove("SimResults")
    }
    i <- i + 1
  }

  if (agree) {
    results <- SimResults$results # nolint: object_name_linter
    plot.sim.seq(results) # nolint: object_name_linter
  } else {
    warning("Results not found. \n")
  }
}


list.sim.seq.scenarios <- function(path) { # nolint: object_name_linter
  # Files in "path" will be scanned and their scenarios listed.

  all.files <- dir(paste(path, "/SimResults", sep = "")) # nolint: object_name_linter
  results <- data.frame(
    filename    = all.files, p = NA, alpha = NA, n0 = NA, n1 = NA,
    corr        = NA,
    u.sd        = NA, u.start.seq = NA, u.ind.fun = "", # nolint: object_name_linter
    distr       = "", descr = ""
  )
  cnt <- 1
  for (fl in all.files) {
    print(fl)
    load(paste(path, "/SimResults/", fl, sep = ""))
    scen <- SimResults$scenario # nolint: object_name_linter
    results$p[cnt]           <- scen$p
    results$alpha[cnt]       <- scen$alpha
    results$n0[cnt]          <- scen$n0
    results$n1[cnt]          <- scen$n1
    results$corr[cnt]        <- ifelse(is.null(scen$corr), 0, scen$corr)
    results$u.sd[cnt]        <- ifelse(is.null(scen$u.sd), NA, scen$u.sd) # nolint: object_name_linter
    results$u.start.seq[cnt] <- ifelse(is.null(scen$u.start.seq), NA, # nolint: object_name_linter
                                       scen$u.start.seq)
    results$u.ind.fun[cnt]   <- ifelse(is.null(scen$u.ind.fun), NA, # nolint: object_name_linter
                                       scen$u.ind.fun)
    results$distr[cnt]       <- scen$distr
    results$descr[cnt]       <- ifelse(is.null(scen$description), NA,
                                       scen$description)
    cnt <- cnt + 1
  }
  results
}
