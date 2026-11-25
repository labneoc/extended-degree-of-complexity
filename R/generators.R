# =============================================================================
# generators.R -- simulated processes used to test the degree of complexity
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
# Requires base R only (tested with R 4.5.3).
# =============================================================================

## Stationary Gaussian process with autocorrelation rfun(tau), tau = 0, 1, ...
## Exact simulation by circulant embedding (Dietrich & Newsam, 1997).
## Returns two independent realisations (real and imaginary parts) and the
## ratio min/max of the embedding eigenvalues (must be >= 0 for exactness).
gauss_circ <- function(N, rfun) {
  r   <- rfun(0:N)
  m   <- 2 * N
  row <- c(r[1:N], r[N + 1], r[N:2])
  lam <- Re(fft(row))
  minlam <- min(lam)
  lam[lam < 0] <- 0
  W <- fft(sqrt(lam / m) * complex(real = rnorm(m), imaginary = rnorm(m)))
  list(x1 = Re(W)[1:N], x2 = Im(W)[1:N], minlam = minlam / max(lam))
}

## Power-law autocorrelation (1 + tau)^(-a); completely monotone, hence
## positive definite.
pl <- function(a) function(t) (1 + t)^(-a)

## Liverani-Saussol-Vaienti intermittent map. Correlations of generic Hoelder
## observables decay as n^(1 - 1/g) (Young 1999; Sarig 2002; Gouezel 2004).
## A uniform jitter in [0, 1e-12] per step keeps the doubling branch from
## collapsing onto dyadic rationals in double precision.
lsv <- function(N, g, burn = 1e4) {
  x   <- runif(1)
  out <- numeric(N)
  eps <- runif(N + burn) * 1e-12
  for (i in 1:(N + burn)) {
    x <- if (x < 0.5) x * (1 + (2 * x)^g) else 2 * x - 1
    x <- x + eps[i]
    if (x >= 1) x <- x - 1 + 1e-9
    if (x <= 0) x <- 1e-9
    if (i > burn) out[i - burn] <- x
  }
  out
}

## Logistic map x -> r x (1 - x).
logistic <- function(N, r, burn = 1e5) {
  x   <- 0.3
  out <- numeric(N)
  for (i in 1:(N + burn)) {
    x <- r * x * (1 - x)
    if (i > burn) out[i - burn] <- x
  }
  out
}

## Fourth-order Runge-Kutta integration of a flow f; step h, one sample every
## `sub` steps; returns the first coordinate.
rk4_flow <- function(N, h, sub, f, s0, burn = 2000) {
  s   <- s0
  out <- numeric(N)
  for (i in 1:(N + burn)) {
    for (j in 1:sub) {
      k1 <- f(s); k2 <- f(s + h / 2 * k1); k3 <- f(s + h / 2 * k2); k4 <- f(s + h * k3)
      s  <- s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    }
    if (i > burn) out[i - burn] <- s[1]
  }
  out
}

## Lorenz (1963) and Roessler (1976) vector fields, classical parameters.
lor <- function(s) c(10 * (s[2] - s[1]), s[1] * (28 - s[3]) - s[2], s[1] * s[2] - 8 / 3 * s[3])
ros <- function(s) c(-s[2] - s[3], s[1] + 0.2 * s[2], 0.2 + s[3] * (s[1] - 5.7))
