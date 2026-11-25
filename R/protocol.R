# =============================================================================
# protocol.R -- finite-data estimation of kappa (original) and kappa* (extended)
#
# Author : Bruno Duarte Gomes (ORCID 0000-0002-3646-7847)
#          Laboratório de Neurofisiologia Eduardo Oswaldo Cruz, Instituto de
#          Ciências Biológicas, Universidade Federal do Pará, 66075-110 Belém,
#          Pará, Brazil.
#          Laboratório de Simulação e Biologia Computacional, Centro de
#          Computação de Alto Desempenho e Inteligência Artificial,
#          Universidade Federal do Pará, 66075-110 Belém, Pará, Brazil.
# Contact: brunodgomes@ufpa.br
#
# Part of: "Extending the degree of complexity: convergence domain, maximal
#           correlation and finite-data estimation"
# Requires base R only; source generators.R and estimators.R first.
# =============================================================================

LAGS   <- unique(round(10^seq(0, log10(20000), by = 0.1)))  # 0.1-decade grid
TAU_LO <- 10                                                # start of fit window
NSUR   <- 40                                                # surrogates per null

## Decision rules (Section 3.4 of the paper), applied in this order:
##  (a) < 5 resolved lags or < 1 decade resolved   -> "unresolved"
##  (b) envelope ratio E(tau_hi)/E(tau_lo) > 0.5    -> "non-decaying"
##  (c) larger R^2 of log-log vs semi-log fit       -> "power" / "exponential"
fit_decay <- function(lags, y, keep) {
  E    <- rev(cummax(rev(abs(y))))          # monotone envelope
  idx  <- which(lags >= TAU_LO)
  run  <- idx[cumprod(keep[idx]) == 1]      # contiguous resolved lags
  span <- if (length(run)) log10(max(lags[run]) / min(lags[run])) else 0
  if (length(run) < 5 || span < 1)
    return(list(cls = "unresolved", alpha = NA, tau_hi = if (length(run)) max(lags[run]) else NA))
  if (E[max(run)] / E[min(run)] > 0.5)
    return(list(cls = "non-decaying", alpha = 0, tau_hi = max(lags[run])))
  lt  <- log(lags[run])
  le  <- log(pmax(E[run], 1e-12))
  r2p <- summary(lm(le ~ lt))$r.squared
  r2e <- summary(lm(le ~ lags[run]))$r.squared
  list(cls    = if (r2e > r2p) "exponential" else "power",
       alpha  = unname(-coef(lm(le ~ lt))[2]),
       tau_hi = max(lags[run]))
}

## kT: original definition (nu >= 0); kX: extended definition (nu > -1).
## Non-decaying cases are reported as non-mixing (kappa undefined); kX = Inf
## is only a flag.
kappas <- function(f) {
  if (f$cls %in% c("unresolved", "exponential")) return(c(kT = 0, kX = 0))
  if (f$cls == "non-decaying") return(c(kT = 1, kX = Inf))
  c(kT = if (f$alpha <= 1) 1 else 1 / f$alpha, kX = 1 / f$alpha)
}

## Full analysis of one series: Pearson + HGR for K in Ks, time-shift null.
analyze <- function(x, label, Ks = c(4, 8, 16)) {
  N  <- length(x)
  xc <- x - mean(x)
  rows <- list(); curves <- list()

  ap  <- acf_fft(x, max(LAGS))[LAGS + 1]
  sp  <- replicate(NSUR, {
    o <- sample(round(N / 4):round(3 * N / 4), 1)
    abs(cor(xc, xc[c((o + 1):N, 1:o)]))
  })
  thr <- mean(sp) + 4 * sd(sp)
  f <- fit_decay(LAGS, ap, abs(ap) > thr); k <- kappas(f)
  rows[[1]]   <- data.frame(label, method = "Pearson", K = NA, cls = f$cls, alpha = f$alpha,
                            tau_hi = f$tau_hi, thr = thr, kappa_T = k["kT"], kappa_ext = k["kX"])
  curves[[1]] <- data.frame(label, method = "Pearson", K = NA, lag = LAGS, val = ap, raw = ap, thr = thr)

  for (K in Ks) {
    h   <- hgr_lags(x, LAGS, K)
    s   <- hgr_null(x, K, NSUR, "shift")
    thr <- mean(s) + 4 * sd(s)
    hd  <- sqrt(pmax(h^2 - mean(s^2), 0))   # debiased rho*_K
    f <- fit_decay(LAGS, hd, h > thr); k <- kappas(f)
    rows[[length(rows) + 1]]     <- data.frame(label, method = "HGR", K = K, cls = f$cls, alpha = f$alpha,
                                               tau_hi = f$tau_hi, thr = thr, kappa_T = k["kT"], kappa_ext = k["kX"])
    curves[[length(curves) + 1]] <- data.frame(label, method = "HGR", K = K, lag = LAGS, val = hd, raw = h, thr = thr)
  }
  list(rows = do.call(rbind, rows), curves = do.call(rbind, curves))
}
