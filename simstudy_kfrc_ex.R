library(MixSim)
library(aricode)
library(parallel)
library(fclust)
library(mvtnorm)
source("simstudy_src.R")

DATA_DIR <- "data"
RESULT_DIR <- "results"
LOG_DIR <- "logs"

for (d in c(DATA_DIR, RESULT_DIR, LOG_DIR)){
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

seed_val <- 1

data_grid <- expand.grid(
  seed = seed_val,
  K = c(2, 4, 6),
  BarOmega = c(0.001, 0.01, 0.05, 0.10),
  sph = c(TRUE, FALSE),
  PiLow = c(1.0, 0.05),
  sn_ratio = c(Inf, 4, 1, 0.25),
  p_total = 10,
  n = 200,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

data_grid$p_noise <- with(data_grid, round(p_total / (1 + sn_ratio)))
data_grid$p_signal <- data_grid$p_total - data_grid$p_noise
data_grid$dataset_id <- sprintf(
  "D-s%02d-K%d-om%.3f-sph%d-pL%.2f-r%s-ps%d-pn%d-n%d",
  data_grid$seed, data_grid$K, data_grid$BarOmega,
  as.integer(data_grid$sph), data_grid$PiLow,
  ifelse(is.finite(data_grid$sn_ratio),
         sprintf("%g", data_grid$sn_ratio), "Inf"),
  data_grid$p_signal, data_grid$p_noise, data_grid$n
)

method_grid <- expand.grid(
  kernel = c("product", "sum"),
  criterion = c("logsum", "sum", "gap"),
  alpha = c(1, 5, 10, 25, Inf),
  fuzz_group = c("half", "quarter", "crt"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)

method_grid$group_id <- with(method_grid, sprintf(
  "G-%s-%s-a%s-%s",
  kernel, criterion,
  ifelse(is.finite(alpha), sprintf("%g", alpha), "Inf"),
  fuzz_group
))

fuzzifier_specs <- list(
  half = list(list(fuzz = "power", m = 2),
              list(fuzz = "exp", m = 2)),
  quarter = list(list(fuzz = "quad", m = 0.5)),
  crt = list(list(fuzz = "crt", m = 2))
)

compute_tau <- function(fuzz_group, K) {
  switch(fuzz_group,
         half = 0.5,
         quarter = 0.25,
         crt = log(K / (K - 1)) / 2,
         stop("Unknown fuzz_group: ", fuzz_group))
}

d_idx <- seq_len(nrow(data_grid))
m_idx <- seq_len(nrow(method_grid))
grid_idx <- expand.grid(d = d_idx, m = m_idx)
run_grid <- cbind(data_grid[grid_idx$d, ],
                  method_grid[grid_idx$m, , drop = FALSE])
rownames(run_grid) <- NULL

cat(sprintf("Datasets:  %d\nGroups:    %d\nUnits:     %d\n",
            nrow(data_grid), nrow(method_grid), nrow(run_grid)))

# Data generation
make_dataset <- function(cfg) {
  set.seed(cfg$seed)
  ms <- MixSim(BarOmega = cfg$BarOmega, K = cfg$K, p = cfg$p_signal,
               sph = as.logical(cfg$sph), hom = TRUE,
               PiLow = cfg$PiLow, resN = 1e5)
  sim <- simdataset(n = cfg$n, Pi = ms$Pi, Mu = ms$Mu, S = ms$S)
  
  X_signal_raw <- sim$X
  X_signal <- scale(X_signal_raw)
  X_noise <- if (cfg$p_noise > 0) {
    matrix(rnorm(cfg$n * cfg$p_noise), nrow = cfg$n, ncol = cfg$p_noise)
  } else {
    matrix(numeric(0), nrow = cfg$n, ncol = 0)
  }
  X <- cbind(X_signal, X_noise)
  colnames(X) <- c(paste0("s", seq_len(cfg$p_signal)),
                   if (cfg$p_noise > 0) paste0("n", seq_len(cfg$p_noise))
                   else character(0))
  list(X = X,
       X_signal_raw = X_signal_raw,
       labels = sim$id,
       signal_idx = seq_len(cfg$p_signal),
       noise_idx = if (cfg$p_noise > 0)
         cfg$p_signal + seq_len(cfg$p_noise) else integer(0),
       mixsim_pars  = list(Pi = ms$Pi,
                           Mu = ms$Mu,
                           S = ms$S,
                           OmegaMap = ms$OmegaMap,
                           BarOmega_realised = ms$BarOmega,
                           BarOmega_requested = cfg$BarOmega),
       cfg = cfg)
}

cache_dataset <- function(cfg) {
  f <- file.path(DATA_DIR, paste0(cfg$dataset_id, ".rds"))
  if (!file.exists(f)) saveRDS(make_dataset(cfg), f)
  f
}

# Diagnostics

importance_weights <- function(bw) { inv <- 1 / bw; inv / sum(inv) }

fuzzy_summaries <- function(U) {
  labels <- max.col(U, ties.method = "first")
  maxU <- apply(U, 1, max)
  ent <- -rowSums(ifelse(U > 0, U * log(U), 0))
  per_cl <- do.call(rbind, lapply(sort(unique(labels)), function(k) {
    idx <- which(labels == k)
    c(cluster = k,
      n = length(idx),
      max_mean = mean(maxU[idx]),
      max_median = median(maxU[idx]),
      max_q25 = unname(quantile(maxU[idx], 0.25)),
      max_q75 = unname(quantile(maxU[idx], 0.75)),
      ent_mean = mean(ent[idx]))
  }))
  list(hard_labels = labels,
       per_cluster = per_cl,
       PC = mean(rowSums(U^2)),
       PE = mean(ent),
       max_global = c(mean = mean(maxU), sd = sd(maxU)))
}

top_eig_summary <- function(eig, K) {
  top <- eig[seq_len(K - 1)]
  gap <- if (length(eig) >= K) eig[K - 1] - eig[K] else NA_real_
  c(sum = sum(top), mean = mean(top), min = min(top),
    max = max(top), gap = gap)
}

# BW diagnostics

build_bw_diag <- function(bw_sel, d, tau, fuzz_group) {
  feasible <- isTRUE(bw_sel$feasible)
  K <- d$cfg$K
  if (!feasible) {
    return(list(
      fuzz_group = fuzz_group,
      tau = tau,
      feasible = FALSE,
      stage1_feasible_count = bw_sel$stage1_feasible_count %||% 0L,
      bw_error = bw_sel$error %||% "infeasible_stage1",
      bw = NA, W = NA, signal_share = NA_real_,
      lambda1 = NA_real_, eigvals_top = NA,
      eigvals_cm1 = setNames(rep(NA_real_, 5),
                             c("sum","mean","min","max","gap")),
      criterion_value = NA_real_, alpha_empirical = NA_real_,
      best_start_id = NA_integer_, iterations = NA_integer_
    ))
  }
  bw <- bw_sel$bw
  W  <- importance_weights(bw)
  best_match <- match(bw_sel$best_start_id,
                      sapply(bw_sel$stage2_results, `[[`, "start_id"))
  list(
    fuzz_group = fuzz_group,
    tau = tau,
    feasible = TRUE,
    stage1_feasible_count = bw_sel$stage1_feasible_count,
    bw_error = NA_character_,
    bw = bw,
    W = W,
    signal_share = sum(W[d$signal_idx]),
    lambda1 = bw_sel$lambda1,
    eigvals_top = bw_sel$eigenvalues[
      seq_len(min(K, length(bw_sel$eigenvalues)))],
    eigvals_cm1 = top_eig_summary(bw_sel$eigenvalues, K),
    criterion_value = bw_sel$criterion_value,
    alpha_empirical = max(bw) / min(bw),
    best_start_id = bw_sel$best_start_id,
    iterations = bw_sel$stage2_results[[best_match]]$iterations
  )
}

# Single unit of work

run_unit <- function(dataset_rds, method) {
  d <- readRDS(dataset_rds)
  X <- d$X
  labels <- d$labels
  K <- d$cfg$K
  actual_sizes <- table(labels)
  true_imbalance <- min(actual_sizes) / max(actual_sizes)
  tau <- compute_tau(method$fuzz_group, K)
  fuzzs <- fuzzifier_specs[[method$fuzz_group]]
  out_files <- vapply(fuzzs, function(fs)
    file.path(RESULT_DIR,
              sprintf("%s__%s__%s.rds",
                      d$cfg$dataset_id, method$group_id, fs$fuzz)),
    character(1))
  if (all(file.exists(out_files))) return(invisible(out_files))
  
  set.seed(d$cfg$seed)
  bw_sel <- tryCatch(
    select_bandwidth(scale(X), c = K, tau = tau, alpha = method$alpha,
                      kernel = method$kernel, criterion = method$criterion,
                      max_iter = 500, lr = 0.2,
                      n_starts = 5, verbose = FALSE),
    error = function(e) list(feasible = FALSE,
                             error = conditionMessage(e))
  )
  bw_diag <- build_bw_diag(bw_sel, d, tau, method$fuzz_group)
  
  P_true <- tryCatch(
    compute_true_posteriors(d$X_signal_raw,
                            d$mixsim_pars$Pi,
                            d$mixsim_pars$Mu,
                            d$mixsim_pars$S),
    error = function(e) NULL
  )
  
  for (i in seq_along(fuzzs)) {
    fs <- fuzzs[[i]]
    res <- c(list(
      dataset_id = d$cfg$dataset_id,
      group_id = method$group_id,
      fuzzifier = fs$fuzz,
      m = fs$m,
      data_cfg = as.list(d$cfg),
      method_cfg = as.list(method),
      timestamp = Sys.time()
    ), bw_diag)
    
    if (!bw_diag$feasible) {
      res$U <- NA
      res$frc_converged <- NA
      res$frc_iter <- NA_integer_
      res$frc_objective <- NA_real_
      res$frc_error <- NA_character_
      res$fuzzy <- NA
      res$ARI <- NA_real_
      res$AMI <- NA_real_
      res$FARI_crisp <- NA_real_
      res$FARI_fuzzy <- NA_real_
      res$true_imbalance <- true_imbalance
      res$est_imbalance <- NA_real_
      saveRDS(res, out_files[i])
      next
    }
    
    set.seed(d$cfg$seed)
    R <- compute_dissimilarity(scale(X), bw_diag$bw, kernel = method$kernel)
    fr <- tryCatch(
      switch(fs$fuzz,
             power = FRC(R^2, k = K, m = fs$m),
             exp   = FRC_exp(R^2, k = K, m = fs$m),
             quad  = FRC_quad(R^2, k = K, m = fs$m),
             crt   = FRC_rc(R^2, k = K, m = fs$m)),
      error = function(e) list(U = NULL, error = conditionMessage(e))
    )
    
    if (is.null(fr$U)) {
      res$U <- NA
      res$frc_converged <- NA
      res$frc_iter <- NA_integer_
      res$frc_objective <- NA_real_
      res$frc_error <- fr$error %||% "unknown"
      res$fuzzy <- NA
      res$ARI <- NA_real_
      res$AMI <- NA_real_
      res$FARI_crisp <- NA_real_
      res$FARI_fuzzy <- NA_real_
      res$true_imbalance <- true_imbalance
      res$est_imbalance <- NA_real_
    } else {
      res$U <- fr$U
      res$frc_converged <- fr$converged
      res$frc_iter <- fr$iter
      res$frc_objective <- fr$objective
      res$frc_error <- NA_character_
      fsum <- fuzzy_summaries(fr$U)
      res$fuzzy <- fsum
      res$ARI <- aricode::ARI(fsum$hard_labels, labels)
      res$AMI <- aricode::AMI(fsum$hard_labels, labels)
      res$FARI_crisp <- tryCatch(
        fclust::ARI.F(labels, fr$U),
        error = function(e) NA_real_
      )
      res$FARI_fuzzy <- if (!is.null(P_true)) {
        tryCatch(fari(P_true, fr$U), error = function(e) NA_real_)
      } else NA_real_
      frc_sizes <- table(fsum$hard_labels)
      res$true_imbalance <- true_imbalance
      res$est_imbalance <- min(frc_sizes) / max(frc_sizes)
    }
    saveRDS(res, out_files[i])
  }
  invisible(out_files)
}

invisible(lapply(split(data_grid, seq_len(nrow(data_grid))), cache_dataset))

method_cols <- names(method_grid)
run_list <- lapply(seq_len(nrow(run_grid)), function(i) {
  list(dataset_rds = file.path(DATA_DIR,
                               paste0(run_grid$dataset_id[i], ".rds")),
       method = as.list(run_grid[i, method_cols]))
})

N_CORES <- min(100L, parallel::detectCores())
cat(sprintf("Running %d units on %d cores\n", length(run_list), N_CORES))

parallel::mclapply(
  run_list,
  function(r) {
    tryCatch(run_unit(r$dataset_rds, r$method),
             error = function(e) {
               cat(sprintf("[%s] %s || %s || %s\n",
                           Sys.time(), basename(r$dataset_rds),
                           r$method$group_id, conditionMessage(e)),
                   file = file.path(LOG_DIR,
                                    sprintf("err_%s.log", Sys.getpid())),
                   append = TRUE)
               NULL
             })
  },
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Aggregation helper
aggregate_results <- function(result_dir = RESULT_DIR, seeds = NULL) {
  files <- list.files(result_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!is.null(seeds)) {
    seed_pat <- paste0("^D-s(",
                       paste(sprintf("%02d", seeds), collapse = "|"),
                       ")-")
    files <- files[grepl(seed_pat, basename(files))]
  }
  rows <- lapply(files, function(f) {
    r <- readRDS(f)
    eig <- r$eigvals_cm1
    if (is.null(eig) || all(is.na(eig)))
      eig <- setNames(rep(NA_real_, 5), c("sum","mean","min","max","gap"))
    data.frame(
      dataset_id = r$dataset_id,
      group_id = r$group_id,
      fuzzifier = r$fuzzifier,
      m = r$m,
      seed = r$data_cfg$seed,
      K = r$data_cfg$K,
      BarOmega = r$data_cfg$BarOmega,
      sph = r$data_cfg$sph,
      PiLow = r$data_cfg$PiLow,
      sn_ratio = r$data_cfg$sn_ratio,
      p_signal = r$data_cfg$p_signal,
      p_noise = r$data_cfg$p_noise,
      kernel = r$method_cfg$kernel,
      criterion = r$method_cfg$criterion,
      alpha = r$method_cfg$alpha,
      fuzz_group = r$method_cfg$fuzz_group,
      tau = r$tau,
      feasible = r$feasible,
      stage1_feasible = r$stage1_feasible_count,
      lambda1 = r$lambda1 %||% NA_real_,
      sum_cm1 = eig["sum"],
      min_cm1 = eig["min"],
      gap_cm1 = eig["gap"],
      crit_value = r$criterion_value %||% NA_real_,
      alpha_emp = r$alpha_empirical %||% NA_real_,
      signal_share = r$signal_share %||% NA_real_,
      ARI = r$ARI %||% NA_real_,
      AMI = r$AMI %||% NA_real_,
      FARI_crisp = r$FARI_crisp %||% NA_real_,
      FARI_fuzzy = r$FARI_fuzzy %||% NA_real_,
      true_imb = r$true_imbalance %||% NA_real_,
      est_imb  = r$est_imbalance %||% NA_real_,
      PC = if (is.list(r$fuzzy)) r$fuzzy$PC else NA_real_,
      PE = if (is.list(r$fuzzy)) r$fuzzy$PE else NA_real_,
      frc_converged   = r$frc_converged %||% NA,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

current_seed <- data_grid$seed[1]
agg <- aggregate_results(seeds = current_seed)
name <- file.path(RESULT_DIR, sprintf("agg_seed%02d.rds", current_seed))
saveRDS(agg, name)