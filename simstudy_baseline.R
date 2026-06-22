library(MixSim)
library(aricode)
library(fclust)
library(mvtnorm)
library(cluster)
library(parallel)
library(dplyr)

source("simstudy_src.R")
source("fuzzy_clust_algs.R")

DATA_DIR <- "data"
RESULT_DIR <- "results"
LOG_DIR <- "logs"
AGG_FILE <- "agg.rds"

OUT <- list(
  fcmmk = file.path(RESULT_DIR, "baseline_fcmmk_res"),
  msfcm = file.path(RESULT_DIR, "baseline_msfcm_res"),
  vfcm = file.path(RESULT_DIR, "baseline_vfcm_res"),
  frc = file.path(RESULT_DIR, "baseline_frc_res"),
  fanny = file.path(RESULT_DIR, "baseline_fanny_res")
)

for (d in c(unlist(OUT), LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

agg <- readRDS(AGG_FILE)
datasets_to_run <- agg %>%
  dplyr::filter(feasible) %>%
  dplyr::pull(dataset_id) %>%
  unique()

cat(sprintf("Datasets to process: %d\n", length(datasets_to_run)))

frc_methods <- list(
  power = list(fn = function(D2, K) FRC(D2, k = K, m = 2), m = 2),
  exp = list(fn = function(D2, K) FRC_exp(D2, k = K, m = 2), m = 2),
  quad = list(fn = function(D2, K) FRC_quad(D2,k = K, m = 0.5), m = 0.5),
  crt = list(fn = function(D2, K) FRC_rc(D2,  k = K, m = 2), m = 2)
)

fanny_methods <- list(
  euclidean = list(metric = "euclidean", m = 2),
  manhattan = list(metric = "manhattan", m = 2)
)

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
  
  # FCM-MK
  out_file <- file.path(OUT$fcmmk, sprintf("%s.rds", dataset_id))
  if (!file.exists(out_file)) {
    set.seed(d$cfg$seed)
    fit <- tryCatch(
      suppressWarnings(fcm_mk(X, c = K, m = 2, max_iter = 100, verbose = FALSE)),
      error = function(e) list(U = NULL, error = conditionMessage(e))
    )
    res <- list(dataset_id = dataset_id, method = "FCM-MK", m = 2,
                data_cfg = as.list(d$cfg), timestamp = Sys.time())
    if (is.null(fit$U)) {
      res$ARI <- res$AMI <- res$FARI_crisp <- res$FARI_fuzzy <- NA_real_
      res$U <- NA; res$converged <- NA; res$iter <- NA_integer_
      res$sigmas <- NA; res$W <- NA
      res$error <- fit$error %||% "unknown"
    } else {
      ev <- evaluate_membership(t(fit$U), labels, P_true)
      res$U <- ev$U; res$ARI <- ev$ARI; res$AMI <- ev$AMI
      res$FARI_crisp <- ev$FARI_crisp; res$FARI_fuzzy <- ev$FARI_fuzzy
      res$converged <- fit$converged; res$iter <- fit$iter
      res$sigmas <- fit$sigmas; res$W <- fit$W; res$J <- fit$J
      res$error <- NA_character_
    }
    saveRDS(res, out_file)
  }
  
  # MSFCM
  out_file <- file.path(OUT$msfcm, sprintf("%s.rds", dataset_id))
  if (!file.exists(out_file)) {
    set.seed(d$cfg$seed)
    fit <- tryCatch(
      suppressWarnings(msfcm(X, c = K, m = 2, max_iter = 200, verbose = FALSE)),
      error = function(e) list(U = NULL, error = conditionMessage(e))
    )
    res <- list(dataset_id = dataset_id, method = "MSFCM", m = 2,
                data_cfg = as.list(d$cfg), timestamp = Sys.time())
    if (is.null(fit$U)) {
      res$ARI <- res$AMI <- res$FARI_crisp <- res$FARI_fuzzy <- NA_real_
      res$U <- NA; res$converged <- NA; res$iter <- NA_integer_
      res$error <- fit$error %||% "unknown"
    } else {
      ev <- evaluate_membership(t(fit$U), labels, P_true)
      res$U <- ev$U; res$ARI <- ev$ARI; res$AMI <- ev$AMI
      res$FARI_crisp <- ev$FARI_crisp; res$FARI_fuzzy <- ev$FARI_fuzzy
      res$converged <- fit$converged; res$iter <- fit$iter
      res$J_fuzzy <- fit$J_fuzzy; res$J_hard <- fit$J_hard
      res$error <- NA_character_
    }
    saveRDS(res, out_file)
  }
  
  # vFCM
  out_file <- file.path(OUT$vfcm, sprintf("%s.rds", dataset_id))
  if (!file.exists(out_file)) {
    set.seed(d$cfg$seed)
    fit <- tryCatch(
      suppressWarnings(vfcm(X, c = K, m0 = 2, a = 0.95, b = 0.05, k = 2,
                            max_iter = 200, verbose = FALSE)),
      error = function(e) list(U = NULL, error = conditionMessage(e))
    )
    res <- list(dataset_id = dataset_id, method = "vFCM",
                m0 = 2, a = 0.95, b = 0.05, k = 2,
                data_cfg = as.list(d$cfg), timestamp = Sys.time())
    if (is.null(fit$U)) {
      res$ARI <- res$AMI <- res$FARI_crisp <- res$FARI_fuzzy <- NA_real_
      res$U <- NA; res$converged <- NA; res$iter <- NA_integer_
      res$m_final <- NA_real_
      res$error <- fit$error %||% "unknown"
    } else {
      ev <- evaluate_membership(t(fit$U), labels, P_true)
      res$U <- ev$U; res$ARI <- ev$ARI; res$AMI <- ev$AMI
      res$FARI_crisp <- ev$FARI_crisp; res$FARI_fuzzy <- ev$FARI_fuzzy
      res$converged <- fit$converged; res$iter <- fit$iter
      res$m_final <- fit$m_final
      res$J_fuzzy <- fit$J_fuzzy; res$J_hard <- fit$J_hard
      res$error <- NA_character_
    }
    saveRDS(res, out_file)
  }
  
  # FRC
  D2 <- as.matrix(dist(X))^2
  for (name in names(frc_methods)) {
    spec <- frc_methods[[name]]
    out_file <- file.path(OUT$frc, sprintf("%s__%s.rds", dataset_id, name))
    if (file.exists(out_file)) next
    set.seed(d$cfg$seed)
    fr <- tryCatch(spec$fn(D2, K),
                   error = function(e) list(U = NULL,
                                            error = conditionMessage(e)))
    res <- list(dataset_id = dataset_id, method = sprintf("FRC_%s", name),
                m = spec$m, data_cfg = as.list(d$cfg), timestamp = Sys.time())
    if (is.null(fr$U)) {
      res$ARI <- res$AMI <- res$FARI_crisp <- res$FARI_fuzzy <- NA_real_
      res$U <- NA; res$converged <- NA
      res$error <- fr$error %||% "unknown"
    } else {
      ev <- evaluate_membership(fr$U, labels, P_true)
      res$U <- ev$U; res$ARI <- ev$ARI; res$AMI <- ev$AMI
      res$FARI_crisp <- ev$FARI_crisp; res$FARI_fuzzy <- ev$FARI_fuzzy
      res$converged <- fr$converged %||% NA
      res$iter <- fr$iter %||% NA_integer_
      res$objective <- fr$objective %||% NA_real_
      res$error <- NA_character_
    }
    saveRDS(res, out_file)
  }
  
  # FANNY
  for (name in names(fanny_methods)) {
    spec <- fanny_methods[[name]]
    out_file <- file.path(OUT$fanny, sprintf("%s__%s.rds", dataset_id, name))
    if (file.exists(out_file)) next
    set.seed(d$cfg$seed)
    fa <- tryCatch(
      cluster::fanny(X, k = K, diss = FALSE, metric = spec$metric,
                     stand = FALSE, memb.exp = spec$m),
      error = function(e) list(membership = NULL,
                               error = conditionMessage(e))
    )
    res <- list(dataset_id = dataset_id, method = sprintf("FANNY_%s", name),
                metric = spec$metric, m = spec$m,
                data_cfg = as.list(d$cfg), timestamp = Sys.time())
    if (is.null(fa$membership)) {
      res$ARI <- res$AMI <- res$FARI_crisp <- res$FARI_fuzzy <- NA_real_
      res$U <- NA
      res$error <- fa$error %||% "unknown"
    } else {
      ev <- evaluate_membership(fa$membership, labels, P_true)
      res$U <- ev$U; res$ARI <- ev$ARI; res$AMI <- ev$AMI
      res$FARI_crisp <- ev$FARI_crisp; res$FARI_fuzzy <- ev$FARI_fuzzy
      res$objective <- fa$objective %||% NA_real_
      res$error <- NA_character_
    }
    saveRDS(res, out_file)
  }
  
  invisible(NULL)
}

n_cores_detected <- parallel::detectCores()
N_CORES <- if (is.na(n_cores_detected)) 1L else min(100L, n_cores_detected)
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

# Aggregate FCM-MK / MSFCM
aggregate_fcm <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(
      dataset_id = r$dataset_id,
      method = r$method,
      m = r$m,
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
      converged = r$converged %||% NA,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

aggregate_vfcm <- function(dir) {
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

aggregate_simple <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    data.frame(
      dataset_id = r$dataset_id,
      method = r$method,
      m = r$m,
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
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

saveRDS(aggregate_fcm(OUT$fcmmk), file.path(RESULT_DIR, "baseline_fcmmk.rds"))
saveRDS(aggregate_fcm(OUT$msfcm), file.path(RESULT_DIR, "baseline_msfcm.rds"))
saveRDS(aggregate_vfcm(OUT$vfcm), file.path(RESULT_DIR, "baseline_vfcm.rds"))
saveRDS(aggregate_simple(OUT$frc), file.path(RESULT_DIR, "baseline_frc.rds"))
saveRDS(aggregate_simple(OUT$fanny), file.path(RESULT_DIR, "baseline_fanny.rds"))