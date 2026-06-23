library(ggplot2)
library(kernlab)
library(RColorBrewer)
library(gridExtra)
source("simstudy_src_upd.R")

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Visualise KFRC output in feature space (any K)
#'
#' For each point, blend the cluster palette colours weighted by U memberships.
#'
#' @param X Original data matrix (post-scaling)
#' @param U Fuzzy membership matrix (n x K)
#' @param bw KFRC-selected bandwidth vector (length p)
#' @param kernel "product" or "sum"
#' @param palette Optional colour vector of length K (defaults to Set1)
#' @param title Optional plot title
visualise_kfrc <- function(X, U, bw, kernel = c("product", "sum"),
                           cluster_palette = NULL, title = NULL) {
  kernel <- match.arg(kernel)
  U <- as.matrix(U)
  n <- nrow(U)
  K_clusters <- ncol(U)
  
  # Build palette
  if (is.null(cluster_palette)) {
    cluster_palette <- if (K_clusters <= 9) {
      RColorBrewer::brewer.pal(max(3, K_clusters), "Set1")[seq_len(K_clusters)]
    } else {
      grDevices::hcl.colors(K_clusters, palette = "Set 2")
    }
  }
  
  palette_rgb <- t(col2rgb(cluster_palette)) / 255
  blended <- U %*% palette_rgb
  blended[blended < 0] <- 0
  blended[blended > 1] <- 1
  point_colours <- rgb(blended[, 1], blended[, 2], blended[, 3])
  
  # Compute kernel matrix using KFRC bandwidths
  K_mat <- compute_kernel_matrix(X, bw, kernel = kernel)
  
  # Centre and normalise (matches your compute_F)
  H <- diag(n) - matrix(1/n, n, n)
  K_centered <- H %*% K_mat %*% H
  K_centered <- (K_centered + t(K_centered)) / 2
  d <- pmax(diag(K_centered), 1e-10)
  D_inv_sqrt <- diag(1 / sqrt(d))
  F_mat <- (1/n) * D_inv_sqrt %*% K_centered %*% D_inv_sqrt
  F_mat <- (F_mat + t(F_mat)) / 2
  
  # Eigendecomposition for projection
  eig <- eigen(F_mat, symmetric = TRUE)
  proj <- eig$vectors[, 1:3] %*% diag(sqrt(pmax(eig$values[1:3], 0)))
  
  # Project to unit sphere
  row_norms <- sqrt(rowSums(proj^2))
  row_norms[row_norms < 1e-10] <- 1
  proj_norm <- proj / row_norms
  
  plot_df <- data.frame(
    PC1 = proj_norm[, 1],
    PC2 = proj_norm[, 2]
  )
  
  # Plot
  p <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
    geom_point(colour = point_colours, size = 2, alpha = 0.85) +
    coord_fixed() +
    annotate("path",
             x = cos(seq(0, 2 * pi, length.out = 200)),
             y = sin(seq(0, 2 * pi, length.out = 200)),
             colour = "grey60", linetype = "dashed") +
    labs(x = expression("kPCA component 1"),
         y = expression("kPCA component 2"),
         title = title %||% sprintf("KFRC fuzzy partition in feature space (K = %d)",
                                    K_clusters)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 9, colour = "grey40"))
  
  return(p)
}

# Configuration

KFRC_DIR <- "results"
DATASETS_DIR <- "real_data"
DATASETS_TO_PLOT <- c("IS", "diabetes", "rice", "seeds", "vehicle", "wine")

# Dataset loader 
load_dataset <- function(name) {
  if (name == "diabetes") {
    d <- readRDS(file.path(DATASETS_DIR, "diabetes.rds"))
    X <- d[, -ncol(d)]; y <- d[, ncol(d)]
  } else if (name == "rice") {
    d <- readRDS(file.path(DATASETS_DIR, "rice.rds"))
    X <- d[, -ncol(d)]; y <- d[, ncol(d)]
  } else if (name == "seeds") {
    data("seeds", package = "datasetsICR")
    X <- seeds[, -ncol(seeds)]; y <- seeds[, ncol(seeds)]
  } else if (name == "wine") {
    data("wine", package = "HDclassif"); X <- wine[, -1]; y <- wine[, 1]
  } else if (name == "vehicle") {
    d <- readRDS(file.path(DATASETS_DIR, "vehicle.rds"))
    X <- d[, -ncol(d)]; y <- d[, ncol(d)]
  } else if (name == "IS") {
    d <- readRDS(file.path(DATASETS_DIR, "IS.rds"))
    X <- d[, -1]; y <- d[, 1]
  } else stop("Unknown dataset: ", name)
  
  X <- as.matrix(X)
  const_cols <- which(apply(X, 2, sd) == 0)
  if (length(const_cols) > 0) {
    cat(sprintf("  Dropping %d constant column(s) for %s\n",
                length(const_cols), name))
    X <- X[, -const_cols, drop = FALSE]
  }
  list(name = name, X = scale(X),
       y = as.integer(as.factor(y)),
       K = length(unique(y)))
}

find_kfrc_file <- function(dataset_name) {
  candidates <- list.files(KFRC_DIR,
                           pattern = sprintf("^kfrc__%s__", dataset_name),
                           full.names = TRUE)
  if (length(candidates) == 0)
    stop(sprintf("No KFRC file for %s found in %s", dataset_name, KFRC_DIR))
  if (length(candidates) > 1)
    warning(sprintf("Multiple KFRC files for %s; using first: %s",
                    dataset_name, basename(candidates[1])))
  candidates[1]
}

# Build plots

plots <- list()
for (ds_name in DATASETS_TO_PLOT) {
  cat(sprintf("\nProcessing %s...\n", ds_name))
  d <- load_dataset(ds_name)
  cat(sprintf("n = %d, p = %d, K = %d\n", nrow(d$X), ncol(d$X), d$K))
  
  kfrc_file <- find_kfrc_file(ds_name)
  r <- readRDS(kfrc_file)
  cat(sprintf("Loaded: %s (FARI = %.3f)\n",
              basename(kfrc_file), r$FARI_crisp))
  
  dataset_label <- if (ds_name == "IS") "Image Segmentation" else tools::toTitleCase(ds_name)
  fuzz_label <- dplyr::recode(r$fuzzifier,
                              power = "Power",
                              exp = "Exponential",
                              quad = "Quadratic",
                              crt = "Complementary Root")
  crit_label <- dplyr::recode(r$criterion,
                              sum = "Sum",
                              logsum = "Log sum",
                              gap = "Eigengap")
  alpha_expr <- if (is.infinite(r$alpha)) bquote(infinity) else r$alpha
  subttl <- bquote(.(fuzz_label) * ", " * .(crit_label) * ", " * alpha == .(alpha_expr))
  
  ttl <- dataset_label
  
  p <- visualise_kfrc(
    X = d$X, U = r$U, bw = r$bw,
    kernel = r$kernel %||% "sum",
    title = ttl
  )
  p <- p + labs(subtitle = subttl)
  plots[[ds_name]] <- p
  cat("Plot for dataset", ds_name, "generated.\n")
}

grid.arrange(grobs = plots, ncol = 3)
