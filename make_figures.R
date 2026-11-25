# =============================================================================
# make_figures.R -- Figures 1-4 of the paper (vector PDF and 600 dpi PNG,
#                   7.2 in wide; Figure 4 is 3.6 in wide)
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
# Usage  : Rscript make_figures.R   (after run_simulation_study.R)
# Needs  : ggplot2 (>= 3.5.0), dplyr, ragg
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("R/generators.R"); source("R/estimators.R"); source("R/protocol.R")
dir.create("figures", showWarnings = FALSE)

FOC <- "#1f5fa8"   # proposed estimator
CMP <- "grey50"    # original definition
DPI <- 600
base_theme <- theme_classic(base_size = 8) +
  theme(plot.tag = element_text(face = "bold", size = 9),
        legend.background = element_blank(), legend.text = element_text(size = 6.5, hjust = 0),
        legend.key.height = unit(8, "pt"), strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 8))
## Every figure is written twice: figures/<name>.pdf (vector, cairo, fonts
## embedded) and figures/<name>.png (600 dpi, ragg).
if (packageVersion("ggplot2") < "3.5.0") base_theme <- base_theme + theme(legend.text.align = 0)
save_fig <- function(name, w, h, draw) {
  grDevices::cairo_pdf(file.path("figures", paste0(name, ".pdf")), width = w, height = h)
  draw(); invisible(dev.off())
  ragg::agg_png(file.path("figures", paste0(name, ".png")), width = w, height = h, units = "in", res = DPI)
  draw(); invisible(dev.off())
}
grid_draw <- function(plots, nr, nc) function() {
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc)))
  for (i in seq_along(plots))
    print(plots[[i]], vp = grid::viewport(layout.pos.row = (i - 1) %/% nc + 1, layout.pos.col = (i - 1) %% nc + 1))
}
## Subscripted symbols are drawn with plotmath so that they do not depend on
## the font having the Unicode subscript digits.
LAB_K8   <- expression("proposed " * kappa[8]^"*" * " (HGR, K = 8)")
LAB_RHO8 <- expression("HGR " * rho[8]^"*" * (tau))
## Legend placed inside the panel; works with ggplot2 3.4 and >= 3.5.
inside <- function(x, y, jx, jy) {
  if (packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = c(jx, jy))
  else theme(legend.position = c(x, y), legend.justification = c(jx, jy))
}
truth_alpha <- function(l) {
  if (grepl("^GaussPL_a", l)) return(as.numeric(sub("GaussPL_a", "", l)))
  if (grepl("^LSV_g", l)) return(1 / as.numeric(sub("LSV_g", "", l)) - 1)
  if (grepl("StochVol|Reparam|Noise", l)) return(as.numeric(sub(".*_a", "", l)))
  NA
}
rows   <- read.csv("results/rows.csv")
curves <- read.csv("results/curves.csv")
rows$alpha_true <- suppressWarnings(sapply(rows$label, truth_alpha))
rows$mk <- ifelse(rows$method == "Pearson", "Pearson", paste0("HGR_K", rows$K))
LEG_T <- "original \u03ba (\u03bd \u2265 0)"; LEG_X <- "extended \u03ba* (\u03bd > \u22121)"

## ---- Figure 1: analytic comparison of the two definitions -------------------
a  <- 10^seq(log10(0.1), log10(4), length.out = 400)
q  <- seq(1.01, 4, length.out = 400)
da <- rbind(data.frame(x = a, k = pmin(1, 1 / a), def = LEG_T), data.frame(x = a, k = 1 / a, def = LEG_X))
dq <- rbind(data.frame(x = q, k = pmin(1, q - 1), def = LEG_T), data.frame(x = q, k = q - 1, def = LEG_X))
cols <- setNames(c(CMP, FOC), c(LEG_T, LEG_X)); lts <- setNames(c("22", "solid"), c(LEG_T, LEG_X))
f1a <- ggplot(da, aes(x, k, colour = def, linetype = def)) +
  annotate("rect", xmin = 0.1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_line(linewidth = 0.7) + scale_colour_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = lts, name = NULL) + scale_x_log10(breaks = c(0.1, 0.3, 1, 3)) +
  coord_cartesian(ylim = c(0, 4)) + labs(x = "tail exponent \u03b1", y = "degree of complexity", tag = "a") +
  base_theme + inside(0.97, 0.97, 1, 1)
f1b <- ggplot(dq, aes(x, k, colour = def, linetype = def)) +
  annotate("rect", xmin = 2, xmax = 4, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_line(linewidth = 0.7) + scale_colour_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = lts, name = NULL) + coord_cartesian(ylim = c(0, 3)) +
  labs(x = "entropic index q", y = "degree of complexity", tag = "b") +
  base_theme + inside(0.03, 0.97, 0, 1)
save_fig("Figure_1", 7.2, 2.8, grid_draw(list(f1a, f1b), 1, 2))

## ---- Figure 2: calibration against known exponents --------------------------
d <- rows %>% filter(grepl("^GaussPL|^LSV", label), mk %in% c("Pearson", "HGR_K8")) %>%
  mutate(system = ifelse(grepl("^Gauss", label), "a   Gaussian processes", "b   LSV intermittent map"),
         estimator = ifelse(mk == "Pearson", LEG_T, "proposed \u03ba*\u2088 (HGR, K = 8)"),
         kappa = ifelse(mk == "Pearson", kappa_T, kappa_ext),
         status = case_when(cls == "power" ~ "power law resolved",
                            cls == "non-decaying" ~ "non-decaying (drawn at 3.9)",
                            TRUE ~ "unresolved or exponential"),
         kappa = ifelse(is.infinite(kappa), 3.9, pmin(kappa, 3.9)))
set.seed(3)
d$xj <- d$alpha_true * exp(runif(nrow(d), -0.04, 0.04)) * ifelse(d$mk == "Pearson", 0.94, 1.06)
ag <- seq(0.25, 2.5, length.out = 400)
tr <- rbind(data.frame(a = ag, k = 1 / ag, which = "true \u03ba* = 1/\u03b1"),
            data.frame(a = ag, k = pmin(1, 1 / ag), which = "true \u03ba = min(1, 1/\u03b1)"))
tr <- merge(tr, data.frame(system = unique(d$system)))
wall <- data.frame(system = c("a   Gaussian processes", "b   LSV intermittent map"), xmin = c(1.05, 1.6), xmax = 2.6)
f2 <- ggplot() +
  geom_rect(data = wall, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey92") +
  geom_line(data = tr, aes(a, k, linetype = which), colour = "black", linewidth = 0.45) +
  geom_point(data = d, aes(xj, kappa, colour = estimator, shape = status), size = 1.9, stroke = 0.7) +
  scale_colour_manual(values = setNames(c(FOC, CMP), c("proposed \u03ba*\u2088 (HGR, K = 8)", LEG_T)),
                      breaks = c(LEG_T, "proposed \u03ba*\u2088 (HGR, K = 8)"), labels = c(LEG_T, LAB_K8), name = NULL) +
  scale_shape_manual(values = c("power law resolved" = 16, "unresolved or exponential" = 4,
                                "non-decaying (drawn at 3.9)" = 2), name = NULL) +
  scale_linetype_manual(values = setNames(c("solid", "22"), c("true \u03ba* = 1/\u03b1", "true \u03ba = min(1, 1/\u03b1)")), name = NULL) +
  scale_x_log10(breaks = c(0.3, 0.5, 1, 2), limits = c(0.24, 2.6)) + coord_cartesian(ylim = c(-0.1, 4)) +
  facet_wrap(~system) + labs(x = "true decay exponent \u03b1", y = "degree of complexity") +
  base_theme + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"),
                     strip.text = element_text(hjust = 0, size = 8.5, face = "bold"))
save_fig("Figure_2", 7.2, 3.8, function() print(f2))

## ---- Figure 3: mechanisms ---------------------------------------------------
cv <- curves %>% group_by(label, method, K) %>% mutate(rep = ceiling(row_number() / length(LAGS))) %>%
  ungroup() %>% filter(rep == 1, method == "Pearson" | K == 8) %>%
  mutate(est = ifelse(method == "Pearson", "Pearson |c(\u03c4)|", "HGR \u03c1*\u2088(\u03c4)"))
pal2 <- c("HGR \u03c1*\u2088(\u03c4)" = FOC, "Pearson |c(\u03c4)|" = CMP)
lx <- scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1k", "10k"))
ly <- scale_y_log10(breaks = c(1e-3, 1e-2, 1e-1, 1), labels = c("0.001", "0.01", "0.1", "1"))
th_in <- base_theme + inside(0.98, 0.98, 1, 1)
panel <- function(lab, tag) {
  dd  <- cv %>% filter(label == lab) %>% mutate(v = pmax(abs(val), 1e-4))
  thr <- dd %>% group_by(est) %>% summarise(thr = mean(thr), .groups = "drop")
  ggplot(dd, aes(lag, v, colour = est)) +
    geom_hline(data = thr, aes(yintercept = thr, colour = est), linetype = "dotted", linewidth = 0.4) +
    geom_line(linewidth = 0.55) + geom_point(size = 0.7) +
    scale_colour_manual(values = pal2, breaks = names(pal2), labels = c(LAB_RHO8, "Pearson |c(\u03c4)|"), name = NULL) +
    lx + ly + coord_cartesian(ylim = c(1e-3, 1)) + labs(x = "lag \u03c4 (samples)", y = "dependence at lag \u03c4", tag = tag) + th_in
}
pa <- panel("StochVol_a0.5", "a")
pc <- panel("Lorenz_x", "c")
db <- cv %>% filter(label %in% c("Reparam_x_a0.6", "Reparam_exp_a0.6")) %>%
  mutate(series = paste0(ifelse(method == "Pearson", "Pearson", "HGR"), ifelse(grepl("exp", label), " on exp(2.5x)", " on x")),
         v = pmax(abs(val), 1e-4))
pb <- ggplot(db, aes(lag, v, colour = series, linetype = series)) + geom_line(linewidth = 0.55) +
  scale_colour_manual(values = c("HGR on x" = FOC, "HGR on exp(2.5x)" = "#7fb0e6", "Pearson on x" = CMP, "Pearson on exp(2.5x)" = "grey20"), name = NULL) +
  scale_linetype_manual(values = c("HGR on x" = "solid", "HGR on exp(2.5x)" = "22", "Pearson on x" = "solid", "Pearson on exp(2.5x)" = "22"), name = NULL) +
  lx + ly + coord_cartesian(ylim = c(1e-3, 1)) + labs(x = "lag \u03c4 (samples)", y = "dependence at lag \u03c4", tag = "b") + th_in
loc_slope <- function(l, v) {
  lv <- log(pmax(v, 1e-12)); ll <- log(l); n <- length(l)
  sapply(1:n, function(i) { j <- max(1, i - 2):min(n, i + 2); if (length(j) < 4) NA else -coef(lm(lv[j] ~ ll[j]))[2] })
}
dd <- cv %>% filter(grepl("^Cutoff", label), method == "HGR") %>% group_by(label) %>%
  mutate(keep = cumprod(raw > thr) == 1, a_loc = loc_slope(lag, val)) %>% filter(keep, lag >= 3) %>% ungroup() %>%
  mutate(L = sub("Cutoff_L", "L = ", label))
tt <- expand.grid(lag = 10^seq(0.5, 4.3, length.out = 200), L = c(200, 2000, 20000)) %>%
  mutate(a_loc = 0.6 * lag / (1 + lag) + lag / L, L = paste0("L = ", L))
pd <- ggplot() + geom_line(data = tt, aes(lag, a_loc, group = L), colour = "black", linewidth = 0.35, linetype = "22") +
  geom_line(data = dd, aes(lag, a_loc, colour = L), linewidth = 0.6) + geom_point(data = dd, aes(lag, a_loc, colour = L), size = 0.8) +
  scale_colour_manual(values = c("L = 200" = "#0b3d73", "L = 2000" = FOC, "L = 20000" = "#8db8e8"), name = NULL) +
  geom_hline(yintercept = 0.6, colour = "grey60", linewidth = 0.3) + lx + coord_cartesian(ylim = c(0, 2.2), xlim = c(3, 2e4)) +
  labs(x = "lag \u03c4 (samples)", y = "local exponent \u03b1(\u03c4)", tag = "d") +
  base_theme + inside(0.02, 0.98, 0, 1)
save_fig("Figure_3", 7.2, 5.6, grid_draw(list(pa, pb, pc, pd), 2, 2))

## ---- Figure 4: null model and the Lorenz artefact ---------------------------
set.seed(2024)
xl   <- rk4_flow(2^18, 0.01, 5, lor, c(1, 1, 20))
hraw <- hgr_lags(xl, LAGS, 8)
s_i  <- hgr_null(xl, 8, NSUR, "iid")
s_s  <- hgr_null(xl, 8, NSUR, "shift")
nulls <- data.frame(null = c("permutation null (i.i.d.)", "time-shift null"),
                    mean = c(mean(s_i), mean(s_s)), thr = c(mean(s_i) + 4 * sd(s_i), mean(s_s) + 4 * sd(s_s)))
write.csv(nulls, "results/lorenz_null_summary.csv", row.names = FALSE)
pal4 <- c("permutation null (i.i.d.)" = "#d95f02", "time-shift null" = FOC)
f4 <- ggplot(data.frame(lag = LAGS, h = hraw), aes(lag, h)) +
  geom_hline(data = nulls, aes(yintercept = thr, colour = null), linewidth = 0.6) +
  geom_line(colour = "black", linewidth = 0.5) + geom_point(size = 0.8) +
  scale_colour_manual(values = pal4, name = "4-s.d. threshold") + lx + ly + coord_cartesian(ylim = c(3e-3, 1)) +
  labs(x = "lag \u03c4 (samples, \u0394t = 0.05)", y = expression(rho[8]^"*" * (tau) * ", Lorenz " * x)) +
  base_theme + inside(0.02, 0.02, 0, 0) + theme(
                     legend.title = element_text(size = 6.5))
save_fig("Figure_4", 3.6, 2.8, function() print(f4))
