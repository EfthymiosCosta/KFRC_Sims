library(aricode)
library(fclust)
library(parallel)
library(dplyr)
library(MixSim)
source("simstudy_src.R")
DATASETS_DIR <- "real_datasets"
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

FUZZIFIERS <- c("power", "exp", "quad", "crt")
CRITERIA <- c("sum", "logsum", "gap")
ALPHAS <- c(1, 5, 10, 25, Inf)

run_kfrc <- function(d, fuzzifier, criterion, alpha) {
  out_file <- file.path(OUT_DIR,
                        sprintf("kfrc__%s__%s__%s__a%s.rds",
                                d$name, fuzzifier, criterion,
                                ifelse(is.infinite(alpha), "Inf",
                                       as.character(alpha))))
  if (file.exists(out_file)) return(invisible(out_file))

  X <- scale(d$X)
  K <- d$ncl
  m <- frc_m_default(fuzzifier)
  tau <- switch(fuzzifier,
                power = 0.5,
                exp = 0.5,
                quad = 0.25,
                crt = log(K / (K - 1)) / 2)

  set.seed(123)
  sel <- tryCatch(
    select_bandwidth2(X, c = K,
                      tau = tau,
                      alpha = alpha,
                      kernel = "sum",
                      criterion = criterion,
                      n_starts = 10,
                      max_iter = 500,
                      verbose = FALSE),
    error = function(e) list(feasible = FALSE,
                             error = conditionMessage(e))
  )

  res <- list(
    dataset = d$name,
    method = "KFRC",
    fuzzifier = fuzzifier,
    criterion = criterion,
    alpha = alpha,
    kernel = "sum",
    tau = tau,
    m = m,
    K = K,
    n = nrow(X),
    p = ncol(X),
    timestamp = Sys.time()
  )

  if (!isTRUE(sel$feasible)) {
    res$feasible <- FALSE
    res$error <- sel$error %||% "Stage 1 infeasible"
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$bw <- NA
    res$lambda1 <- NA_real_
    res$U <- NA
    saveRDS(res, out_file)
    return(invisible(out_file))
  }

  res$feasible <- TRUE
  res$bw <- sel$bw
  res$lambda1 <- sel$lambda1
  res$eigenvalues <- sel$eigenvalues
  res$criterion_value <- sel$criterion_value

  D2 <- compute_dissimilarity_sq(X, sel$bw, kernel = "sum")
  fit <- tryCatch(
    suppressWarnings(frc_call(fuzzifier, D2, K)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )

  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$uniformity <- NA_real_
    res$U <- NA
    res$frc_converged <- NA
    res$frc_error <- fit$error %||% "frc failed"
  } else {
    ev <- evaluate(fit$U, d$y)
    res$U <- fit$U
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$uniformity <- ev$uniformity
    res$frc_converged <- fit$converged
    res$frc_iter <- fit$iter
    res$frc_error <- NA_character_
  }
  saveRDS(res, out_file)
  invisible(out_file)
}

process_dataset_kfrc <- function(name) {
  cat(sprintf("\n=== %s (KFRC) ===\n", name))
  d <- tryCatch(load_dataset(name), error = function(e) {
    cat(sprintf("  Failed to load: %s\n", conditionMessage(e)))
    NULL
  })
  if (is.null(d)) return(invisible(NULL))

  cat(sprintf("n = %d, p = %d, K = %d\n",
              nrow(d$X), ncol(d$X), d$ncl))

  for (fz in FUZZIFIERS)
    for (cr in CRITERIA)
      for (al in ALPHAS) {
        cat(sprintf("KFRC | %s | %s | a=%s\n", fz, cr,
                    ifelse(is.infinite(al), "Inf", as.character(al))))
        run_kfrc(d, fz, cr, al)
      }
  invisible(NULL)
}

# Run in parallel
N_CORES <- min(100L, parallel::detectCores())
cat(sprintf("Running %d datasets across %d cores\n",
            length(DATASET_NAMES), N_CORES))

parallel::mclapply(
  DATASET_NAMES,
  function(name) {
    tryCatch(process_dataset_kfrc(name),
             error = function(e) {
               cat(sprintf("[%s] %s || %s\n",
                           Sys.time(), name, conditionMessage(e)),
                   file = file.path(LOG_DIR,
                                    sprintf("err_kfrc_%s.log", Sys.getpid())),
                   append = TRUE)
               NULL
             })
  },
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Aggregate
files <- list.files(OUT_DIR, pattern = "^kfrc__", full.names = TRUE)
agg_kfrc <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  data.frame(
    dataset = r$dataset,
    method = "KFRC",
    fuzzifier = r$fuzzifier,
    criterion = r$criterion,
    alpha = r$alpha,
    kernel = r$kernel,
    K = r$K,
    n = r$n,
    p = r$p,
    ARI = r$ARI %||% NA_real_,
    AMI = r$AMI %||% NA_real_,
    FARI_crisp = r$FARI_crisp %||% NA_real_,
    uniformity = r$uniformity %||% NA_real_,
    lambda1 = r$lambda1 %||% NA_real_,
    feasible = r$feasible %||% NA,
    stringsAsFactors = FALSE
  )
}))
saveRDS(agg_kfrc, "results/realdata_kfrc_agg.rds")