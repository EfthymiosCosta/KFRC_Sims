library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

AGG_DIR <- "results"
agg_files <- list.files(AGG_DIR, pattern = "^agg_seed\\d+\\.rds$", full.names = TRUE)
agg <- do.call(rbind, lapply(agg_files, readRDS))
saveRDS(agg, file = "results/agg.rds")

# Helpers
fmt_snr <- function(x) {
  lab <- ifelse(is.infinite(x), "SNR = Inf", sprintf("SNR = %g", x))
  factor(lab, levels = c("SNR = 0.25", "SNR = 1", "SNR = 4", "SNR = Inf"))
}
fmt_c <- function(x) factor(sprintf("c = %d", x), levels = sprintf("c = %d", sort(unique(x))))

my_theme <- theme_bw() +
  theme(plot.title = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = "grey60"),
        strip.text = element_text(colour = "grey20", face = "bold"))

# Oracle picks
oracle_picks <- agg %>%
  filter(feasible) %>%
  group_by(dataset_id) %>%
  slice_max(FARI_fuzzy, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(alpha_lab = factor(ifelse(is.infinite(alpha), "Inf",
                                   as.character(alpha)),
                            levels = c("1", "5", "10", "25", "Inf")),
         snr_lab = fmt_snr(sn_ratio),
         c_lab   = fmt_c(K),
         crit_lab = factor(criterion, levels = c("logsum", "sum", "gap")))

# Panel A: Oracle alpha
pA <- ggplot(oracle_picks, aes(x = alpha_lab, fill = alpha_lab)) +
  geom_bar() +
  facet_grid(c_lab ~ snr_lab) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, guide = "none") +
  labs(x = expression("Oracle " * alpha), y = "Count") +
  my_theme

# Panel B: Oracle criterion
pB <- ggplot(oracle_picks, aes(x = crit_lab, fill = crit_lab)) +
  geom_bar() +
  facet_grid(c_lab ~ snr_lab) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Oracle criterion", y = "Count") +
  my_theme

# Panel C: W_eff boxplots
weff_df <- readRDS("results/weff_df.rds")
pC <- weff_df %>%
  mutate(alpha_lab = factor(ifelse(is.infinite(alpha), "Inf",
                                   as.character(alpha)),
                            levels = c("1", "5", "10", "25", "Inf")),
         snr_lab = fmt_snr(sn_ratio),
         c_lab   = fmt_c(K),
         crit_lab = factor(criterion, levels = c("sum", "logsum", "gap"))) %>%
  ggplot(aes(x = alpha_lab, y = W_eff, fill = crit_lab)) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "grey60") +
  geom_hline(aes(yintercept = p_signal), linetype = "dashed", colour = "darkgreen") +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  facet_grid(c_lab ~ snr_lab) +
  scale_y_log10() +
  scale_fill_manual(values = c(sum = "#F8766D", logsum = "#619CFF", gap = "#00BA38"),
                    name = "Criterion") +
  labs(x = expression(alpha * " bound"),
       y = expression(W[eff] * " (log scale)")) +
  my_theme +
  theme(legend.position = "bottom")

# Combine into one figure with patchwork labels
final_plot <- (pA + pB) / pC +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")

# NEW
baseline_fanny <- readRDS('results/baseline_fanny.rds')
baseline_fcmmk <- readRDS('results/baseline_fcmmk.rds')
baseline_frc <- readRDS('results/baseline_frc.rds')
baseline_msfcm <- readRDS('results/baseline_msfcm.rds')
baseline_vfcm <- readRDS('results/baseline_vfcm.rds')
# 
# agg$U <- (1 - agg$PC) / (1 - 1 / agg$K)
# feas <- aggregate(feasible ~ fuzzifier + alpha, data = agg, FUN = mean)
# cfg <- aggregate(cbind(ARI, AMI, FARI_crisp, FARI_fuzzy, U) ~ sn_ratio + K + fuzzifier + kernel + criterion + alpha,
#                  data = agg, FUN = mean, na.rm = TRUE, na.action = na.pass)
# cfg$fuzzifier <- factor(cfg$fuzzifier,
#                         levels = c("power", "exp", "quad", "crt"),
#                         labels = c("Power", "Exponential", "Quadratic", "Comp. root"))
# 
# cols <- c("method", "sn_ratio", "K", "ARI", "AMI", "FARI_crisp", "FARI_fuzzy", "uniformity")
# base <- rbind(baseline_fanny[cols], baseline_fcmmk[cols], baseline_frc[cols], baseline_msfcm[cols], baseline_vfcm[cols])
# names(base)[names(base) == "uniformity"] <- "U"
# bm <- aggregate(cbind(ARI, AMI, FARI_crisp, FARI_fuzzy, U) ~ method + sn_ratio + K,
#                 data = base, FUN = mean, na.rm = TRUE, na.action = na.pass)
# 
# frc <- bm[grepl("^FRC_", bm$method), ]
# frc$fuzzifier <- factor(sub("FRC_", "", frc$method),
#                         levels = c("power", "exp", "quad", "crt"),
#                         labels = c("Power", "Exponential", "Quadratic", "Comp. root"))
# oth <- bm[!grepl("^FRC_", bm$method), ]
# 
# box_metric <- function(metric) {
#   ggplot(cfg, aes(x = fuzzifier, y = .data[[metric]])) +
#     geom_boxplot(outlier.size = 0.5) +
#     geom_hline(data = oth, aes(yintercept = .data[[metric]], colour = method), linetype = "dashed") +
#     geom_point(data = frc, aes(shape = "FRC"), colour = "red", size = 2.5) +
#     scale_shape_manual(values = 4) +
#     facet_grid(K ~ sn_ratio,
#                labeller = labeller(K = function(x) paste("c =", x),
#                                    sn_ratio = function(x) paste("SNR =", x))) +
#     labs(x = NULL, y = metric, colour = NULL, shape = NULL) +
#     theme_bw() +
#     theme(legend.position = "bottom")
# }
# 
# box_metric("FARI_fuzzy")

# Sensitivity of mean FARI_fuzzy against alpha per fuzzifier, SNR, and number of clusters
agg$U <- (1 - agg$PC) / (1 - 1 / agg$K)
ds <- aggregate(cbind(ARI, AMI, FARI_crisp, FARI_fuzzy, U) ~ dataset_id + sn_ratio + K + fuzzifier + alpha,
                data = agg, FUN = mean, na.rm = TRUE, na.action = na.pass)
plot_alpha <- function(metric) {
  m <- aggregate(ds[[metric]], by = ds[c("sn_ratio", "K", "fuzzifier", "alpha")],
                 FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                     se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                                     n = sum(!is.na(x))))
  m <- do.call(data.frame, m)
  m <- m[m$x.n >= 10, ]
  m$alpha <- factor(m$alpha, levels = c(1, 5, 10, 25, Inf), labels = c("1", "5", "10", "25", "Inf"))
  m$fuzzifier <- factor(m$fuzzifier,
                        levels = c("power", "exp", "quad", "crt"),
                        labels = c("Power", "Exponential", "Quadratic", "Comp. Root"))
  ylab <- switch(metric, FARI_fuzzy = expression(FARI[fuzzy]), FARI_crisp = expression(FARI[crisp]), metric)
  pd <- position_dodge(width = 0.3)
  ggplot(m, aes(x = alpha, y = x.mean, colour = fuzzifier, group = fuzzifier)) +
    geom_line(position = pd) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = x.mean - 1.96 * x.se, ymax = x.mean + 1.96 * x.se),
                  width = 0.2, position = pd) +
    facet_grid(K ~ sn_ratio,
               labeller = labeller(K = function(x) paste("c =", x),
                                   sn_ratio = function(x) paste("SNR =", x))) +
    labs(x = expression(alpha), y = ylab, colour = NULL) +
    scale_x_discrete(labels = c("1", "5", "10", "25", expression(infinity))) +
    theme_bw() +
    theme(legend.position = "bottom")
}

plot_alpha("FARI_fuzzy")

# ablation study - repeated-measures ANOVA
key <- c("dataset_id", "sn_ratio", "K", "fuzzifier", "kernel", "criterion")
fa <- agg[agg$feasible, ]
fa <- fa[order(fa$alpha), ]
rule <- fa[!duplicated(fa[key]), ]

make_long <- function(metric) {
  kf <- aggregate(rule[[metric]], by = rule[c("dataset_id", "sn_ratio", "K", "fuzzifier")], FUN = mean, na.rm = TRUE)
  names(kf)[ncol(kf)] <- "y"
  kf$method <- paste0("KFRC_", kf$fuzzifier)
  cols <- c("dataset_id", "sn_ratio", "K", "method", metric)
  base <- rbind(baseline_fanny[cols], baseline_fcmmk[cols], baseline_frc[cols], baseline_msfcm[cols], baseline_vfcm[cols])
  names(base)[ncol(base)] <- "y"
  long <- rbind(kf[c("dataset_id", "sn_ratio", "K", "method", "y")], base)
  long <- long[is.finite(long$y), ]
  cnt <- table(long$dataset_id)
  long[long$dataset_id %in% names(cnt)[cnt == length(unique(long$method))], ]
}

rm_anova <- function(metric) {
  long <- make_long(metric)
  abl <- long[grepl("^(KFRC|FRC)_", long$method), ]
  abl$fuzzifier <- factor(sub(".*_", "", abl$method), levels = c("power", "exp", "quad", "crt"))
  abl$method <- factor(sub("_.*", "", abl$method), levels = c("FRC", "KFRC"))
  abl$dataset_id <- factor(abl$dataset_id)
  cells <- split(abl, list(abl$sn_ratio, abl$K), drop = TRUE)
  res <- do.call(rbind, lapply(cells, function(d) {
    d$dataset_id <- droplevels(d$dataset_id)
    s <- summary(aov(y ~ method * fuzzifier + Error(dataset_id / (method * fuzzifier)), data = d))
    s <- s[names(s) != "Error: dataset_id"]
    tabs <- do.call(rbind, lapply(s, function(st) {
      tab <- st[[1]]
      rn <- trimws(rownames(tab))
      ss_res <- tab[rn == "Residuals", "Sum Sq"]
      eff <- tab[rn != "Residuals", ]
      data.frame(effect = trimws(rownames(eff)),
                 F = eff[["F value"]],
                 p = eff[["Pr(>F)"]],
                 eta2p = eff[["Sum Sq"]] / (eff[["Sum Sq"]] + ss_res))
    }))
    data.frame(sn_ratio = d$sn_ratio[1], K = d$K[1], n = nlevels(d$dataset_id), tabs)
  }))
  rownames(res) <- NULL
  m <- aggregate(y ~ sn_ratio + K + method + fuzzifier, data = abl, FUN = mean)
  m$combo <- paste(m$method, m$fuzzifier, sep = "_")
  m <- reshape(m[c("sn_ratio", "K", "combo", "y")], idvar = c("sn_ratio", "K"), timevar = "combo", direction = "wide")
  m[-(1:2)] <- round(m[-(1:2)], 3)
  list(anova = res, means = m)
}

abl <- rm_anova("FARI_fuzzy")
abl$anova
abl$means

# Interaction plots comparing performance of FRC and KFRC
plot_ablation <- function(metric) {
  long <- make_long(metric)
  abl <- long[grepl("^(KFRC|FRC)_", long$method), ]
  abl$fuzzifier <- factor(sub(".*_", "", abl$method),
                          levels = c("power", "exp", "quad", "crt"),
                          labels = c("Power", "Exponential", "Quadratic", "Comp. root"))
  abl$method <- factor(sub("_.*", "", abl$method), levels = c("FRC", "KFRC"))
  m <- aggregate(abl$y, by = abl[c("sn_ratio", "K", "method", "fuzzifier")],
                 FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  m <- do.call(data.frame, m)
  ylab <- switch(metric, FARI_fuzzy = expression(FARI[fuzzy]), FARI_crisp = expression(FARI[crisp]), metric)
  pd <- position_dodge(width = 0.2)
  ggplot(m, aes(x = fuzzifier, y = x.mean, colour = method, group = method)) +
    geom_line(position = pd) +
    geom_point(position = pd) +
    geom_errorbar(aes(ymin = x.mean - 1.96 * x.se, ymax = x.mean + 1.96 * x.se),
                  width = 0.15, position = pd) +
    facet_grid(K ~ sn_ratio,
               labeller = labeller(K = function(x) paste("c =", x),
                                   sn_ratio = function(x) paste("SNR =", x))) +
    labs(x = NULL, y = ylab, colour = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_ablation("FARI_fuzzy")