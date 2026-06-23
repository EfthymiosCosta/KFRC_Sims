library(aricode)
library(fclust)
library(cluster)
library(parallel)
library(dplyr)
library(MixSim)

source("simstudy_src.R")
source("fuzzy_clust_algs.R")

DATASETS_DIR <- "real_data"
OUT_DIR <- "results"
LOG_DIR <- "realdata_logs"
for (d in c(OUT_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

load_dataset <- function(name) {
  if (name == "diabetes") {
    d <- readRDS(file.path(DATASETS_DIR, "diabetes.rds"))
    X <- d[, -ncol(d)]
    y <- d[, ncol(d)]
  } else if (name == "rice") {
    d <- readRDS(file.path(DATASETS_DIR, "rice.rds"))
    X <- d[, -ncol(d)]
    y <- d[, ncol(d)]
  } else if (name == "seeds") {
    if (!requireNamespace("datasetsICR", quietly = TRUE))
      stop("Package 'datasetsICR' required for the seeds dataset.")
    data("seeds", package = "datasetsICR")
    X <- seeds[, -ncol(seeds)]
    y <- seeds[, ncol(seeds)]
  } else if (name == "wine") {
    if (!requireNamespace("HDclassif", quietly = TRUE))
      stop("Package 'HDclassif' required for the wine dataset.")
    data("wine", package = "HDclassif")
    X <- wine[, -1]
    y <- wine[, 1]
  } else if (name == "vehicle") {
    d <- readRDS(file.path(DATASETS_DIR, "vehicle.rds"))
    X <- d[, -ncol(d)]; y <- d[, ncol(d)]
  } else if (name == "IS") {
    d <- readRDS(file.path(DATASETS_DIR, "IS.rds"))
    X <- d[, -1]
    y <- d[, 1]
  } else {
    stop(sprintf("Unknown dataset: %s", name))
  }
  
  X <- as.matrix(X)
  const_cols <- which(apply(X, 2, sd) == 0)
  if (length(const_cols) > 0) {
    cat(sprintf("[%s] Dropping %d constant column(s): %s\n",
                name, length(const_cols),
                paste(colnames(X)[const_cols], collapse = ", ")))
    X <- X[, -const_cols, drop = FALSE]
  }
  
  list(name = name,
       X = X,
       y = as.integer(as.factor(y)),
       ncl = length(unique(y)))
}

DATASET_NAMES <- c("diabetes", "rice", "seeds",
                   "wine", "vehicle", "IS")

uniformity_score <- function(U, K) {
  PC <- mean(rowSums(U^2))
  (1 - PC) / (1 - 1 / K)
}

leading_eig_euclid <- function(X) {
  R <- as.matrix(dist(X))
  D2 <- R^2
  n <- nrow(D2)
  H <- diag(n) - matrix(1/n, n, n)
  B <- -0.5 * H %*% D2 %*% H
  d <- diag(B)
  D_inv_sqrt <- diag(1 / sqrt(d))
  F_mat <- (1/n) * D_inv_sqrt %*% B %*% D_inv_sqrt
  if (requireNamespace("RSpectra", quietly = TRUE))
    return(RSpectra::eigs_sym(F_mat, k = 1)$values[1])
  return(eigen(F_mat, symmetric = TRUE, only.values = TRUE)$values[1])
}

evaluate <- function(U, y) {
  hard <- max.col(U, ties.method = "first")
  list(
    ARI = aricode::ARI(hard, y),
    AMI = aricode::AMI(hard, y),
    FARI_crisp = tryCatch(fclust::ARI.F(y, U),
                          error = function(e) NA_real_),
    uniformity = uniformity_score(U, ncol(U))
  )
}

frc_call <- function(fuzzifier, D2, K) {
  switch(fuzzifier,
         power = FRC(D2, k = K, m = 2),
         exp = FRC_exp(D2, k = K, m = 2),
         quad = FRC_quad(D2, k = K, m = 0.5),
         crt = FRC_rc(D2, k = K, m = 2),
         stop("Unknown fuzzifier"))
}

frc_m_default <- function(fuzzifier) {
  switch(fuzzifier, power = 2, exp = 2, quad = 0.5, crt = 2)
}

# Methods
run_frc_uniform <- function(d, fuzzifier) {
  out_file <- file.path(OUT_DIR,
                        sprintf("frc__%s__%s.rds", d$name, fuzzifier))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl
  D2 <- as.matrix(dist(X))^2

  set.seed(123)
  fit <- tryCatch(
    suppressWarnings(frc_call(fuzzifier, D2, K)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )

  res <- list(
    dataset = d$name, method = "FRC",
    fuzzifier = fuzzifier, m = frc_m_default(fuzzifier),
    K = K, n = nrow(X), p = ncol(X),
    timestamp = Sys.time()
  )
  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$error <- fit$error %||% "frc failed"
  } else {
    ev <- evaluate(fit$U, d$y)
    res$U <- fit$U
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$converged <- fit$converged
    res$iter <- fit$iter
    res$error <- NA_character_
  }
  saveRDS(res, out_file); invisible(out_file)
}

run_fanny <- function(d, metric) {
  out_file <- file.path(OUT_DIR,
                        sprintf("fanny__%s__%s.rds", d$name, metric))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl

  set.seed(123)
  fit <- tryCatch(
    suppressWarnings(cluster::fanny(X, k = K, diss = FALSE,
                                    metric = metric, stand = FALSE,
                                    memb.exp = 2)),
    error = function(e) list(membership = NULL,
                             error = conditionMessage(e))
  )

  res <- list(dataset = d$name, method = "FANNY",
              metric = metric, m = 2, K = K,
              n = nrow(X), p = ncol(X),
              timestamp = Sys.time())
  if (is.null(fit$membership)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$error <- fit$error %||% "fanny failed"
  } else {
    ev <- evaluate(fit$membership, d$y)
    res$U <- fit$membership
    res$ARI <- ev$ARI; res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$error <- NA_character_
  }
  saveRDS(res, out_file); invisible(out_file)
}

run_msfcm <- function(d) {
  out_file <- file.path(OUT_DIR, sprintf("msfcm__%s.rds", d$name))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl

  set.seed(123)
  fit <- tryCatch(
    suppressWarnings(msfcm(X, c = K, m = 2, max_iter = 200, verbose = FALSE)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )

  res <- list(dataset = d$name, method = "MSFCM",
              m = 2, K = K, n = nrow(X), p = ncol(X),
              timestamp = Sys.time())
  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$error <- fit$error %||% "msfcm failed"
  } else {
    U_nc <- t(fit$U)
    ev <- evaluate(U_nc, d$y)
    res$U <- U_nc
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$converged <- fit$converged
    res$iter <- fit$iter
    res$error <- NA_character_
  }
  saveRDS(res, out_file); invisible(out_file)
}

run_vfcm <- function(d) {
  out_file <- file.path(OUT_DIR, sprintf("vfcm__%s.rds", d$name))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl

  set.seed(123)
  fit <- tryCatch(
    suppressWarnings(vfcm(X, c = K, m0 = 2, a = 0.95, b = 0.05, k = 2,
                          max_iter = 200, verbose = FALSE)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )
  res <- list(dataset = d$name, method = "vFCM",
              m0 = 2, K = K, n = nrow(X), p = ncol(X),
              timestamp = Sys.time())
  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$error <- fit$error %||% "vfcm failed"
  } else {
    U_nc <- t(fit$U)
    ev <- evaluate(U_nc, d$y)
    res$U <- U_nc
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$m_final <- fit$m_final
    res$converged <- fit$converged
    res$iter <- fit$iter
    res$error <- NA_character_
  }
  saveRDS(res, out_file); invisible(out_file)
}

run_fcmmk <- function(d) {
  out_file <- file.path(OUT_DIR, sprintf("fcmmk__%s.rds", d$name))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl

  set.seed(123)
  fit <- tryCatch(
    suppressWarnings(fcm_mk(X, c = K, m = 2, max_iter = 100, verbose = FALSE)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )

  res <- list(dataset = d$name, method = "FCM-MK",
              m = 2, K = K, n = nrow(X), p = ncol(X),
              timestamp = Sys.time())
  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$error <- fit$error %||% "fcm-mk failed"
  } else {
    U_nc <- t(fit$U)
    ev <- evaluate(U_nc, d$y)
    res$U <- U_nc
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$sigmas <- fit$sigmas
    res$W <- fit$W
    res$converged <- fit$converged
    res$iter <- fit$iter
    res$error <- NA_character_
  }
  saveRDS(res, out_file)
  invisible(out_file)
}

FUZZIFIERS <- c("power", "exp", "quad", "crt")

process_dataset_baselines <- function(name) {
  cat(sprintf("\n=== %s (baselines) ===\n", name))
  d <- tryCatch(load_dataset(name), error = function(e) {
    cat(sprintf("Failed to load: %s\n", conditionMessage(e)))
    NULL
  })
  if (is.null(d)) return(invisible(NULL))

  X_scaled <- scale(d$X)
  l1_eu <- tryCatch(leading_eig_euclid(X_scaled),
                    error = function(e) NA_real_)
  saveRDS(list(dataset = d$name, n = nrow(X_scaled),
               p = ncol(X_scaled), K = d$ncl,
               leading_eig_euclid = l1_eu),
          file.path(OUT_DIR, sprintf("descript__%s.rds", d$name)))

  cat(sprintf("n = %d, p = %d, K = %d, lambda_1(Euclid) = %.4f\n",
              nrow(X_scaled), ncol(X_scaled), d$ncl, l1_eu))

  for (fz in FUZZIFIERS) {
    cat(sprintf("FRC | %s\n", fz))
    run_frc_uniform(d, fz)
  }
  for (mt in c("euclidean", "manhattan")) {
    cat(sprintf("FANNY | %s\n", mt))
    run_fanny(d, mt)
  }
  cat("MSFCM\n")
  run_msfcm(d)
  cat("vFCM\n")
  run_vfcm(d)
  cat("FCM-MK\n")
  run_fcmmk(d)

  invisible(NULL)
}

# Run in parallel

N_CORES <- min(100L, parallel::detectCores())
cat(sprintf("Running %d datasets across %d cores\n",
            length(DATASET_NAMES), N_CORES))

parallel::mclapply(
  DATASET_NAMES,
  function(name) {
    tryCatch(process_dataset_baselines(name),
             error = function(e) {
               cat(sprintf("[%s] %s || %s\n",
                           Sys.time(), name, conditionMessage(e)),
                   file = file.path(LOG_DIR,
                                    sprintf("err_baselines_%s.log",
                                            Sys.getpid())),
                   append = TRUE)
               NULL
             })
  },
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Aggregation

aggregate_baselines <- function() {
  files <- list.files(OUT_DIR, pattern = "^(frc|fanny|msfcm|vfcm|fcmmk)__",
                      full.names = TRUE)
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(
      dataset = r$dataset,
      method = r$method,
      fuzzifier = r$fuzzifier %||% NA_character_,
      metric = r$metric %||% NA_character_,
      m = r$m %||% NA_real_,
      K = r$K %||% NA_integer_,
      n = r$n %||% NA_integer_,
      p = r$p %||% NA_integer_,
      ARI = r$ARI %||% NA_real_,
      AMI = r$AMI %||% NA_real_,
      FARI_crisp = r$FARI_crisp %||% NA_real_,
      uniformity = r$uniformity %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

agg_baselines <- aggregate_baselines()
saveRDS(agg_baselines, "results/realdata_baselines_agg.rds")