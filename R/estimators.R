# =============================================================================
# estimators.R -- Pearson autocorrelation and binned maximal correlation
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

## Sample autocorrelation up to maxlag by FFT, unbiased normalisation 1/(N - tau).
acf_fft <- function(x, maxlag) {
  N <- length(x)
  x <- x - mean(x)
  f <- fft(c(x, numeric(N)))
  a <- Re(fft(Mod(f)^2, inverse = TRUE))[1:(maxlag + 1)] / (2 * N)
  a <- a / (N - 0:maxlag)
  a / a[1]
}

## Rank-based assignment to K equiprobable cells (quantile partition).
binify <- function(x, K) {
  as.integer(ceiling(rank(x, ties.method = "first") / length(x) * K))
}

## Binned Hirschfeld-Gebelein-Renyi maximal correlation between two cell
## sequences: second singular value of Q_ij = P_ij / sqrt(p_i p_j).
hgr1 <- function(b1, b2, K) {
  P <- matrix(tabulate((b1 - 1L) * K + b2, K * K), K, K, byrow = TRUE) / length(b1)
  Q <- P / sqrt(outer(rowSums(P), colSums(P)))
  svd(Q, nu = 0, nv = 0)$d[2]
}

## rho*_K(tau) on a set of lags.
hgr_lags <- function(x, lags, K) {
  b <- binify(x, K)
  N <- length(b)
  sapply(lags, function(t) hgr1(b[1:(N - t)], b[(1 + t):N], K))
}

## Null distributions for the HGR statistic.
## iid  : random permutation of the partner sequence (destroys autocorrelation;
##        shown in the paper to produce false long-range dependence).
## shift: circular shift by a random offset in [N/4, 3N/4] (keeps each
##        series' own autocorrelation; used in the final protocol).
hgr_null <- function(x, K, nsur = 40, type = c("shift", "iid")) {
  type <- match.arg(type)
  b <- binify(x, K)
  N <- length(b)
  replicate(nsur, {
    if (type == "shift") {
      o <- sample(round(N / 4):round(3 * N / 4), 1)
      hgr1(b, b[c((o + 1):N, 1:o)], K)
    } else {
      hgr1(b[1:(N - 1)], b[2:N][sample(N - 1)], K)
    }
  })
}
