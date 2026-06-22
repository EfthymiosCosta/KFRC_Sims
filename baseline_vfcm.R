library(aricode)
library(fclust)
library(mvtnorm)
library(parallel)
library(dplyr)
source("simstudy_src_upd.R")
source("vFCM.R")

DATA_DIR <- "sim_bw_selection/data"
OUT_DIR <- "sim_bw_selection/baseline_vfcm_res"
LOG_DIR <- "sim_bw_selection/baseline_vfcm_logs"
for (d in c(OUT_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

agg <- readRDS("agg.rds")
datasets_to_run <- agg %>%
  filter(feasible) %>%
  pull(dataset_id) %>%
  unique()

cat(sprintf("Datasets to process: %d\n", length(datasets_to_run)))

evaluate_membership <- function(U, labels, P_true) {
  hard <- max.col(U, ties.method = "first")
  list(
    U = U,
    hard = hard,
    ARI = aricode::ARI(hard, labels),
    AMI = aricode::AMI(hard, labels),
    FARI_crisp = tryCatch(fclust::ARI.F(labels, U),
                          error = function(e) NA_real_),
    FARI_fuzzy = if (!is.null(P_true))
      tryCatch(fari(P_true, U), error = function(e) NA_real_)
    else NA_real_
  )
}

run_one_dataset <- function(dataset_id) {
  out_file <- file.path(OUT_DIR, sprintf("%s.rds", dataset_id))
  if (file.exists(out_file)) return(invisible(out_file))

  d <- tryCatch(
    readRDS(file.path(DATA_DIR, paste0(dataset_id, ".rds"))),
    error = function(e) NULL
  )
  if (is.null(d)) return(invisible(NULL))

  X <- scale(d$X)
  labels <- d$labels
  K <- d$cfg$K

  P_true <- tryCatch(
    compute_true_posteriors(d$X_signal_raw,
                            d$mixsim_pars$Pi,
                            d$mixsim_pars$Mu,
                            d$mixsim_pars$S),
    error = function(e) NULL
  )

  set.seed(d$cfg$seed)
  fit <- tryCatch(
    suppressWarnings(vfcm(X, c = K, m0 = 2, a = 0.95, b = 0.05, k = 2,
                          max_iter = 200, verbose = FALSE)),
    error = function(e) list(U = NULL, error = conditionMessage(e))
  )

  res <- list(
    dataset_id = dataset_id,
    method = "vFCM",
    m0 = 2,
    a = 0.95,
    b = 0.05,
    k = 2,
    data_cfg = as.list(d$cfg),
    timestamp = Sys.time()
  )

  if (is.null(fit$U)) {
    res$ARI <- NA_real_
    res$AMI <- NA_real_
    res$FARI_crisp <- NA_real_
    res$FARI_fuzzy <- NA_real_
    res$U <- NA
    res$converged <- NA
    res$iter <- NA_integer_
    res$m_final <- NA_real_
    res$error <- fit$error %||% "unknown"
  } else {
    U_nc <- t(fit$U)
    ev <- evaluate_membership(U_nc, labels, P_true)
    res$U <- ev$U
    res$ARI <- ev$ARI
    res$AMI <- ev$AMI
    res$FARI_crisp <- ev$FARI_crisp
    res$FARI_fuzzy <- ev$FARI_fuzzy
    res$converged <- fit$converged
    res$iter <- fit$iter
    res$m_final <- fit$m_final
    res$J_fuzzy <- fit$J_fuzzy
    res$J_hard <- fit$J_hard
    res$error <- NA_character_
  }

  saveRDS(res, out_file)
  invisible(out_file)
}

N_CORES <- min(100L, parallel::detectCores())
cat(sprintf("Running on %d cores\n", N_CORES))

parallel::mclapply(
  datasets_to_run,
  function(d_id) {
    tryCatch(run_one_dataset(d_id),
             error = function(e) {
               cat(sprintf("[%s] %s || %s\n",
                           Sys.time(), d_id, conditionMessage(e)),
                   file = file.path(LOG_DIR,
                                    sprintf("err_%s.log", Sys.getpid())),
                   append = TRUE)
               NULL
             })
  },
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

aggregate_baseline <- function(dir = OUT_DIR) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(
      dataset_id = r$dataset_id,
      method = r$method,
      m0 = r$m0,
      seed = r$data_cfg$seed,
      K = r$data_cfg$K,
      BarOmega = r$data_cfg$BarOmega,
      sph = r$data_cfg$sph,
      PiLow = r$data_cfg$PiLow,
      sn_ratio = r$data_cfg$sn_ratio,
      p_signal = r$data_cfg$p_signal,
      p_noise = r$data_cfg$p_noise,
      ARI = r$ARI %||% NA_real_,
      AMI = r$AMI %||% NA_real_,
      FARI_crisp = r$FARI_crisp %||% NA_real_,
      FARI_fuzzy = r$FARI_fuzzy %||% NA_real_,
      iter = r$iter %||% NA_integer_,
      m_final = r$m_final %||% NA_real_,
      converged = r$converged %||% NA,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

agg_vfcm <- aggregate_baseline()
saveRDS(agg_vfcm, "baseline_vfcm.rds")
cat(sprintf("Saved %d rows to baseline_vfcm.rds\n", nrow(agg_vfcm)))
