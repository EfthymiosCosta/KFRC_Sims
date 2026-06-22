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