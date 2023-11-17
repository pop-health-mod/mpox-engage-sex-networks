# Mpox SEIR model ----
mpox_mod <- function(init_pop = init_pop,
                     contact = contact,
                     beta = 0.7, alpha = 1 / 5.1, gamma = 1 / 5, omega = 0.5,
                     report_frac = 0.7, report_delay = 1 / 2, diag_ass = TRUE,
                     start = 0, end = 90, dt = 0.01,
                     calibration = FALSE) {
  
  # model time step  
  niter <- (end - start) / dt + 1
  time <- seq(start, end, dt)
  K <- ncol(init_pop)
  
  # SEIR model and indices
  X <- array(data = 0, dim = c(niter, 6, K))
  s <- 1
  e1 <- 2
  e2 <- 3
  i1 <- 4
  i2 <- 5
  r <- 6
  
  # accounting of reported cases
  R <- array(data = 0, dim = c(niter, 2, K))
  X[1, ,] <- init_pop
  
  # proportional mixing
  mix_prp <- matrix(data = NA, nrow = K, ncol = K)
  for (k in 1:K) {
    for (k_ in 1:K) {
      mix_prp[k, k_] <- (contact[k_] * sum(X[1,, k_])) / sum(contact * colSums(X[1,,]))
    }
  }
  
  # assortative mixing
  mix_ass <- diag(K)
  if(diag_ass == FALSE) {
    if (!((K / 5)%%1 == 0)) { print("error, number of groups need to be a factor of 5 with diag_ass = FALSE") }
    mix_ass <- matrix(0, nrow = K, ncol = K)
    index <- data.frame(lw = seq(1, K, 5),
                        up = seq(5, K, 5))
    for (i in 1:nrow(index)) {
      mix_ass[index$lw[i]:index$up[i], index$lw[i]:index$up[i]] <- 1 / 5
    }
  }
  
  # final mixing
  mix <- omega * mix_ass + (1 - omega) * mix_prp
  
  # Euler loop
  cases <- rep(0, niter)
  inc <- rep(0, niter)
  Rt <- rep(0, niter)
  if (calibration == FALSE) {
    ngm <- sweep(mix, MARGIN = 1, contact, "*")
    Rt[1] <- eigen(ngm * beta * (1 / gamma))$values[1]
    # Rt[1] <- ngm_erlang2_fun(beta = beta, alpha = alpha, gamma = gamma,  K = K, mix = mix)    
  }
  
  for (t in 2:niter) {
    
    # group-specific prevalence
    prv <- matrix(data = 0, nrow = 1, ncol = K)
    for (k in 1:K) {   
      prv[k] <- sum(X[t - 1, c(i1, i2), k]) / sum(X[t - 1, , k])
    } 
    
    # foi
    trans_rate <- matrix(data = 0, nrow = K, ncol = K)
    for (k in 1:K) {
      for (k_ in 1:K) {
        trans_rate[k, k_] <- contact[k] * mix[k, k_] * prv[k_] * beta
      }
    }
    foi_t <- rowSums(trans_rate)
    
    # looping over the risk groups (ODE)
    for (k in 1:K) {
      # Susceptible
      X[t, s, k] <-  X[t - 1, s, k] + dt * (- foi_t[k] * X[t - 1, s, k])   
      # Exposed
      X[t, e1, k] <-  X[t - 1, e1, k] + dt * (+ foi_t[k] * X[t - 1, s, k]
                                              - alpha * 2 * X[t - 1, e1, k])   
      X[t, e2, k] <-  X[t - 1, e2, k] + dt * (+ alpha * 2 * X[t - 1, e1, k]
                                              - alpha * 2 * X[t - 1, e2, k])        
      # Infectious
      X[t, i1, k] <-  X[t - 1, i1, k] + dt * (+ alpha * 2 * X[t - 1, e2, k]
                                              - gamma * 2 * X[t - 1, i1, k])  
      X[t, i2, k] <-  X[t - 1, i2, k] + dt * (+ gamma * 2 * X[t - 1, i1, k]
                                              - gamma * 2 * X[t - 1, i2, k])  
      # Recovered
      X[t, r, k] <-  X[t - 1, r, k] + dt * (+ gamma * 2 * X[t - 1, i2, k])        
      
      # Cases that become infectious/symptomatic
      R[t, 1, k] <-  R[t - 1, 1, k] + dt * (+ alpha * 2 * X[t - 1, e2, k]
                                            - report_delay * R[t - 1, 1, k]) 
      # Cases that are reported
      R[t, 2, k] <-  R[t - 1, 2, k] + dt * (+ report_frac * report_delay * R[t - 1, 1, k]) 
    } 
    if (any(X[t,,] < 0)) { 
      print("negative compartments, parameter set rejected")
      return(list(cases = rep(-Inf, niter), time = time, X = X, R = R, mix = mix, mix_prp = mix_prp, Rt = Rt))
      break }
    inc[t] <- sum(foi_t * X[t - 1, s, ])
    cases[t] <- sum(report_frac * report_delay * R[t - 1, 1, ]) # per day
    
    if (calibration == FALSE) {    
      Rt[t] <- eigen(sweep(ngm , MARGIN = 1, X[t - 1, s, ] / colSums(X[t - 1, ,]), "*") * beta * (1 / gamma))$values[1]
      #Rt[t] <- ngm_erlang2_fun(beta = beta, alpha = alpha, gamma = gamma,
      #                           K = K, mix = mix,
      #                          prp_sus = X[t - 1, s, ] / colSums(X[t - 1, ,]))
    }
  }
  
  return(list(cases = cases, time = time, X = X, R = R, mix = mix, mix_prp = mix_prp, Rt = Rt, inc = inc))
}

# ---- log-likelihood ----
## TODO CLEAN UP FUNCTIONS
## TODO ADD COMMENTS TO DOCUMENT
prior_dens <- function(theta) {
  # priors
  # plogis(qlogis(0.5) + c(-1, 1) * qnorm(0.95) * 2) # beta (transmission parameter)
  # plogis(qlogis(0.5) + c(-1, 1) * qnorm(0.95) * 1) # omega (mixing)
  # plogis(qlogis(0.8) + c(-1, 1) * qnorm(0.95) * 1) # reporting fraction
  # 1 / exp(log(1/3) + c(1, -1) * qnorm(0.95) * 0.5) # duration of infectiousness
  log_prior <- dnorm(theta[1], mean = qlogis(0.5), sd = 2, log = TRUE) +
               dnorm(theta[2], mean = qlogis(0.5), sd = 1, log = TRUE) +
               dnorm(theta[3], mean = qlogis(0.8), sd = 1, log = TRUE) +
               dnorm(theta[4], mean = log(1 / 3), sd = 0.5, log = TRUE) 
  return(log_prior)
  
}

llk <- function(theta, fixed_par, likdat, calibration = TRUE) {
  # prior theta
  log_prior <- prior_dens(theta)
  
  # getting predictions
  val <- mpox_mod(init_pop = fixed_par$init_pop,
                  contact = fixed_par$contact,
                  beta = plogis(theta[1]), alpha = fixed_par$alpha, gamma = exp(theta[4]), omega = plogis(theta[2]),
                  report_frac = plogis(theta[3]), report_delay = fixed_par$report_delay,
                  start = fixed_par$start, end = ceiling(max(likdat$time_intro)), dt = fixed_par$dt,
                  calibration = TRUE)
  mod_inc <- val$cases[val$time %in% likdat$time_intro]
  
  # likelihood
  log_lik <- sum(dpois(x = likdat$incidence, lambda = mod_inc, log = TRUE))
  
  # posterior
  post_llk <- log_prior + log_lik
  return(post_llk)
}

# get cri function
getci <- function(df) {
  n_day <- ncol(df)
  # Function for the confidence interval of the estimates
  ci <- apply(X = df, MARGIN = 2,
              FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  lower <- ci[1, ]
  upper <- ci[2, ]
  df1 <- data.frame(lower, upper)
  return(df1)
}

# simulation to obtain CrI
simul_fun <- function(hessian, theta, fixed_par, sim = 1000, 
                      SIR = TRUE, nsir = 1000, likdat = NULL,
                      with_replacement = TRUE, track_sims = FALSE,
                      parallel = TRUE) {
  # From the hessian, simulate the model
  vcova <- Matrix::solve(-hessian)
  
  # Test if positive semi-definite
  eS <- eigen(vcova, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -1e-06 * abs(ev[1L]))) {
    vcova <- Matrix::nearPD(vcova, corr = FALSE)$mat
    vcova <- matrix(vcova@x, 
                    nrow = vcova@Dim[1], 
                    ncol = vcova@Dim[2] )
  }
  
  # Laplace approximation only
  if (SIR == FALSE) {
    samp <- mvtnorm::rmvnorm(n = sim, theta, vcova)
  }
  
  if (SIR == TRUE) {
    log_sum_exp <- function(x) {
      xmax <- which.max(x)
      log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax] }
    
    if (is.null(likdat)) { print('SIR needs to include <likdat>'); break }
    # Sample parameters from multivariate t-distribution to get thicker tails
    par_sir <- mvtnorm::rmvt(n = nsir, delta = theta, sigma = as.matrix(vcova), df = 2)
    # We add the mode to the resampled values (to be sure that it is included in the CI)
    par_sir <- rbind(par_sir, theta)
    prp_dens <- mvtnorm::dmvt(par_sir, theta, as.matrix(vcova), df = 2, log = TRUE)
    
    if (parallel == FALSE)  {
      # You might get some warnings, that's OK, the function will return NA
      # and we replace them with -Inf in the llk. We suppress them.
      loglikelihood_sir <- suppressWarnings(
        apply(
          par_sir, MARGIN = 1, FUN = function(x) {
            llk(x, fixed_par, likdat, calibration = TRUE)
          }
        )
      )
    } else {
      require(doParallel)
      cl <- parallel::makeCluster(parallel::detectCores() - 1)
      doParallel::registerDoParallel(cl)
      parallel::clusterExport(cl, c("fixed_par", "likdat", "llk", "prior_dens", "mpox_mod"))
      # You might get some warnings, that's OK, the function will return NA
      # and we replace them with -Inf in the llk. We suppress them.
      loglikelihood_sir <- foreach::foreach(i = 1:nrow(par_sir), .combine = "c") %dopar% {
        par_sir_i <- par_sir[i, ]
        llk(par_sir_i, fixed_par, likdat, calibration = TRUE)  }
      parallel::stopCluster(cl)
    }
    
    loglikelihood_sir[!is.finite(loglikelihood_sir)] <- -Inf
    dens_ratio <- loglikelihood_sir - prp_dens
    wgt <- exp(dens_ratio - log_sum_exp(dens_ratio))
    resampleid <- sample(nrow(par_sir), size = sim, replace = with_replacement, prob = wgt)
    samp <- par_sir[resampleid,, drop = FALSE]
    nunique <- length(unique(resampleid))
    max_wgt <- max(wgt, na.rm = TRUE)
    if (max_wgt > 0.1) { print(paste('Caution, maximum weight is ', round(max_wgt, 2), " consider increasing sim or sampling without replacement (latter option is not the best)")) }
    print(paste(nunique, 'unique parameter sets resampled'))
  }     
  
  val <- matrix(data = NA, nrow = sim, ncol = (fixed_par$end - fixed_par$start) / fixed_par$dt + 1)
  val_r <- val
  for (s in 1:sim) {
    if(track_sims & s %% 100 == 0){ print(sprintf("%s out of %s sims", s, sim)) }
    tmp <- mpox_mod(init_pop = fixed_par$init_pop,
                    contact = fixed_par$contact,
                    beta = plogis(samp[s, 1]), alpha = fixed_par$alpha, gamma = exp(samp[s, 4]), omega = plogis(samp[s, 2]),
                    report_frac = plogis(samp[s, 3]), report_delay = fixed_par$report_delay,
                    start = fixed_par$start, end = fixed_par$end, dt = fixed_par$dt)
    val[s, ] <- tmp$cases
    val_r[s, ] <- tmp$Rt
  }
  
  # format results
  res <- getci(val)
  res_rt <- getci(val_r)
  result <- data.frame(time = seq(fixed_par$start, fixed_par$end, fixed_par$dt),
                       lci = res$lower,
                       uci = res$upper,
                       rt_lci = res_rt$lower,
                       rt_uci = res_rt$upper)
  par_ci <- getci(samp)
  posterior_par <- data.frame(names = c("transmission parameter (beta)", 
                                        "assortativity (omega)", 
                                        "reporting fraction", 
                                        "duration infectiousness (1/gamma)"),
                              lci = c(plogis(par_ci$lower[1:3]), 1 / exp(par_ci$upper[4])),
                              uci = c(plogis(par_ci$upper[1:3]), 1 / exp(par_ci$lower[4])))
  
  return(list(result = result, posterior_ci = posterior_par, posterior_samples = samp))
}


# NGM computation of R_0 / R_eff ----
# simplified approach, equivalent to deriving the transmission and transition matrices
r0_fun <- function(fit, theta, fixed_par, post_samples,
                   grp_immunity = NULL, compute_ci = TRUE) {
  
  ## full NGM to get R0
  if (length(grp_immunity) == 0) {
    mat <- sweep(fit$mix, MARGIN = 1, fixed_par$contact, "*")
    ngm <- mat * plogis(theta[1]) * (1 / exp(theta[4]))
    
    # point estimate
    r0_pt <- eigen(ngm)$values[1]
    
    # ci
    if (compute_ci) {
      r0 <- rep(0, nrow(post_samples))
      for (i in 1:nrow(post_samples)) {
        mix_ass <- diag(length(fit$X[1, 1, ]))
        mix <- plogis(post_samples[i, 2]) * mix_ass + (1 - plogis(post_samples[i, 2])) * fit$mix_prp
        
        mat <- sweep(mix, MARGIN = 1, fixed_par$contact, "*")
        ngm <- mat * plogis(post_samples[i, 1]) * (1 / exp(post_samples[i, 4]))
        r0[i] <- eigen(ngm)$value[1]
      }
    }
    
  ## assume certain proportion immune and get R_eff
  } else {
    # proportion immune
    n_grp <- length(fit$X[1, 1, ])
    if (grp_immunity > n_grp) { print("error, too many groups"); break }
    
    # set their rates to 0 (immune so do not contribute to transmission)
    contact_0 <- fixed_par$contact
    contact_0[(n_grp - grp_immunity + 1):n_grp] <- 0
    
    mat <- sweep(fit$mix, MARGIN = 1, contact_0, "*")
    ngm <- mat * plogis(theta[1]) * (1 / exp(theta[4]))
    
    # point estimate
    r0_pt <- eigen(ngm)$values[1]
    
    # ci
    if (compute_ci) {
      r0 <- rep(0, nrow(post_samples))
      for (i in 1:nrow(post_samples)) {
        mix_ass <- diag(length(fit$X[1, 1, ]))
        mix <- plogis(post_samples[i, 2]) * mix_ass + (1 - plogis(post_samples[i, 2])) * fit$mix_prp
        
        mat <- sweep(mix, MARGIN = 1, contact_0, "*")
        ngm <- mat * plogis(post_samples[i, 1]) * (1 / exp(post_samples[i, 4]))
        r0[i] <- eigen(ngm)$value[1]
      }
    }
  }
  
  # if not computing the CIs, use the point estimate as placeholder
  if (compute_ci) {
    r0_ci <- quantile(r0, probs = c(0.025, 0.975))
    val <- data.frame(pt_est = r0_pt, lci = r0_ci[1], uci = r0_ci[2])
  } else {
    val <- data.frame(pt_est = r0_pt, lci = r0_pt, uci = r0_pt)
    r0 <- r0_pt
  }
  
  return(val)
}

# Diekman ways of doing things
# The 4*4 T matrix adapted to the Erlang distrubiton
ngm_erlang2_fun <- function(beta = 0.9, alpha = 1/5, gamma = 1/5, 
                            K = NULL, mix = mix, 
                            prp_sus = NULL) {
  Tr <- matrix(data = 0, ncol = 4, nrow = 4)
  Sg <- matrix(data = 0, ncol = 4, nrow = 4)
  Sg[1, 1] <- -2 * alpha
  Sg[2, 1] <- +2 * alpha
  Sg[2, 2] <- -2 * alpha
  Sg[3, 2] <- +2 * alpha
  Sg[3, 3] <- -2 * gamma
  Sg[4, 3] <- +2 * gamma
  Sg[4, 4] <- -2 * gamma
  
  Sig <- matrix(data = 0, ncol = K * nrow(Sg), nrow = K * nrow(Sg))
  index <- data.frame(lw = seq(1, ncol(Sig), ncol(Sg)),
                      hg = seq(ncol(Sg), ncol(Sig) + 1, ncol(Sg)))
  for (i in 1:nrow(index)) {
    Sig[index$lw[i]:index$hg[i], index$lw[i]:index$hg[i]] <- Sg
  }
  
  Tra <- matrix(data = 0, ncol = K * nrow(Tr), nrow = K * nrow(Tr))
  index <- data.frame(lw = seq(1, ncol(Tra), ncol(Tr)),
                      hg = seq(ncol(Tr), ncol(Tra) + 1, ncol(Tr)))
  if (length(prp_sus) == 0) {
    for (i in 1:K) {
      for (j in 1:K) {
        infections <- contact[i] * mix[i, j] * beta
        tmp <- matrix(data = 0, nrow = 4, ncol = 4)
        tmp[1, 3:4] <- rep(infections, 2)
        Tra[index$lw[i]:index$hg[i], index$lw[j]:index$hg[j]] <- tmp
      }
    }        
  } else {
    for (i in 1:K) {
      for (j in 1:K) {
        infections <- prp_sus[i] * contact[i] * mix[i, j] * beta
        tmp <- matrix(data = 0, nrow = 4, ncol = 4)
        tmp[1, 3:4] <- rep(infections, 2)
        Tra[index$lw[i]:index$hg[i], index$lw[j]:index$hg[j]] <- tmp
      }
    }
  }
  
  invSig <- solve(Sig)
  r0 <- eigen((-Tra) %*% invSig)$values[1]
  return(r0)
}
