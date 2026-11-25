# =============================================================================
# run_simulation_study.R -- generates every series and estimates kappa, kappa*
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
# Usage  : Rscript run_simulation_study.R   (from the repository root)
# Output : results/rows.csv   (one row per series x estimator)
#          results/curves.csv (dependence curves on the lag grid)
# Runtime: 7 to 13 min on a single core.
# =============================================================================

source("R/generators.R")
source("R/estimators.R")
source("R/protocol.R")
dir.create("results", showWarnings = FALSE)

set.seed(777)
N   <- 2^20
ALL <- list()
add <- function(r) ALL[[length(ALL) + 1]] <<- r

## Gaussian processes with r(tau) = (1 + tau)^(-alpha); 4 series per alpha
for (a in c(0.3, 0.6, 0.9, 1.2, 1.6, 2.0)) for (k in 1:2) {
  g <- gauss_circ(N, pl(a))
  add(analyze(g$x1, sprintf("GaussPL_a%.1f", a)))
  add(analyze(g$x2, sprintf("GaussPL_a%.1f", a)))
}
## Power law with exponential cut-off
for (L in c(200, 2000, 20000)) {
  g <- gauss_circ(N, function(t) (1 + t)^(-0.6) * exp(-t / L))
  add(analyze(g$x1, sprintf("Cutoff_L%d", L)))
  add(analyze(g$x2, sprintf("Cutoff_L%d", L)))
}
## Stochastic volatility: zero linear autocorrelation, long memory in |x|
g <- gauss_circ(N, pl(0.5))
for (z in list(g$x1, g$x2)) add(analyze(exp(0.75 * z) * rnorm(N), "StochVol_a0.5"))
## Monotone relabelling y = exp(2.5 x)
g <- gauss_circ(N, pl(0.6))
add(analyze(g$x1, "Reparam_x_a0.6")); add(analyze(exp(2.5 * g$x1), "Reparam_exp_a0.6"))
add(analyze(g$x2, "Reparam_x_a0.6")); add(analyze(exp(2.5 * g$x2), "Reparam_exp_a0.6"))
## Additive white noise, SNR = 1
g <- gauss_circ(N, pl(0.6))
for (z in list(g$x1, g$x2)) add(analyze(z + rnorm(N), "Noise_SNR1_a0.6"))
## Short-memory controls
for (k in 1:2) {
  add(analyze(as.numeric(arima.sim(list(ar = 0.9), N)), "AR1_phi0.9"))
  add(analyze(rnorm(N), "WhiteNoise"))
}
## LSV intermittent map, alpha_true = 1/g - 1
for (gm in c(0.3, 0.45, 0.6, 0.75)) for (k in 1:3) add(analyze(lsv(N, gm), sprintf("LSV_g%.2f", gm)))
## Logistic map: fully chaotic and Feigenbaum accumulation point
for (k in 1:2) add(analyze(logistic(N, 4) + runif(N) * 1e-12, "Logistic_r4"))
add(analyze(logistic(N, 3.5699456718709449), "Logistic_Feigenbaum"))
## Flows (N = 2^18). Lorenz sampled every 0.05 time units, Roessler every 0.1
NF <- 2^18
for (k in 1:2) {
  add(analyze(rk4_flow(NF, 0.01, 5, lor, c(1, 1, 20) + rnorm(3)), "Lorenz_x"))
  add(analyze(rk4_flow(NF, 0.02, 5, ros, c(1, 1, 0) + rnorm(3) * 0.1), "Rossler_x"))
}

write.csv(do.call(rbind, lapply(ALL, `[[`, "rows")),   "results/rows.csv",   row.names = FALSE)
write.csv(do.call(rbind, lapply(ALL, `[[`, "curves")), "results/curves.csv", row.names = FALSE)
