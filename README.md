# Extending the degree of complexity

Code for the paper *Extending the degree of complexity: convergence domain, maximal correlation and finite-data estimation*, submitted to *Chaos, Solitons & Fractals*.

**Author:** Bruno Duarte Gomes (ORCID [0000-0002-3646-7847](https://orcid.org/0000-0002-3646-7847)), brunodgomes@ufpa.br
Laboratório de Neurofisiologia Eduardo Oswaldo Cruz, Instituto de Ciências Biológicas, Universidade Federal do Pará, 66075-110 Belém, Pará, Brazil.
Laboratório de Simulação e Biologia Computacional, Centro de Computação de Alto Desempenho e Inteligência Artificial, Universidade Federal do Pará, 66075-110 Belém, Pará, Brazil.

## What it does

Tsallis (Particles 9 (2026) 33) defined the degree of complexity as κ = 1/(ν_max + 1), where ν_max is the moment order at which ∫ ξ^ν c(ξ) dξ diverges. For simulated processes with known decay of correlations, this code computes:

- the original κ from the Pearson autocorrelation (ν ≥ 0);
- the extended κ* = 1/α (ν > −1);
- κ*_K from the binned Hirschfeld–Gebelein–Rényi maximal correlation ρ*_K(τ), K ∈ {4, 8, 16}.

## Layout

| File | Content |
|---|---|
| `R/generators.R` | Gaussian processes by circulant embedding, LSV intermittent map, logistic map, Lorenz and Rössler flows (RK4) |
| `R/estimators.R` | FFT autocorrelation, quantile binning, maximal correlation ρ*_K, permutation and time-shift nulls |
| `R/protocol.R` | Lag grid, decision rules (unresolved / non-decaying / power / exponential), κ and κ* |
| `run_simulation_study.R` | Generates every series in the paper and writes `results/rows.csv` and `results/curves.csv` |
| `make_figures.R` | Figures 1–4 as vector PDF and 600 dpi PNG in `figures/`; also writes `results/lorenz_null_summary.csv` |
| `check_reproduction.R` | Compares a new run with the stored output, series by series |
| `results/` | Output of the run used in the paper (Tables 2 and 3, Figures 2 and 3) |
| `figures/` | The figures as submitted |

## Running

Requires R ≥ 4.0. The simulation uses base R only. The figures need `ggplot2` (3.4 or later), `dplyr` and `ragg`.

```
Rscript run_simulation_study.R   # 7 to 13 min on a single core
Rscript make_figures.R
```

Run both from the repository root. The seeds are fixed: `set.seed(777)` for the study and `set.seed(2024)` for the Lorenz series in Figure 4. `run_simulation_study.R` overwrites `results/`; keep a copy of the stored files if you want to compare.

## Reproducibility

The stored output in `results/` was produced with R 4.5.3. A full re-run with R 4.3.3 on Linux (x86-64) gave the following, as reported by `check_reproduction.R`:

- every Gaussian, cut-off, stochastic-volatility, relabelling, additive-noise, AR(1), white-noise, logistic, Lorenz and Rössler series was identical to the stored output, row by row (196 of 196 rows);
- `results/lorenz_null_summary.csv` (Figure 4) was identical;
- the LSV intermittent map at γ = 0.30 was identical, while 35 of the 36 rows for γ = 0.45, 0.60 and 0.75 differed.

The LSV map is chaotic. The most likely cause of the difference is the term (2x)^γ, which is evaluated by the C library's `pow` function; its result can differ in the last bit between platforms and R builds, and the chaotic dynamics amplifies such a difference over 2^20 iterations. The LSV trajectories therefore depend on the machine, although their statistical properties do not. The random-number stream is unaffected, which is why every series generated after the LSV map is again identical. The LSV values in the paper are those in `results/rows.csv`. With only three series per γ, a re-run on another platform gives different individual LSV estimates; the paper treats this replication level as a limitation.

## License

MIT, see `LICENSE`.

## Citation

If you use this code, please cite the paper (reference to be added on publication) and Tsallis, C. (2026). Complexity: what is it? *Particles*, 9, 33. https://doi.org/10.3390/particles9020033
