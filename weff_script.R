suppressPackageStartupMessages({
  library(parallel)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

extract_one <- function(f) {
  r <- readRDS(f)
  if (!isTRUE(r$feasible) || is.null(r$W) || all(is.na(r$W))) return(NULL)
  data.frame(
    seed = r$data_cfg$seed,
    K = r$data_cfg$K,
    BarOmega = r$data_cfg$BarOmega,
    sn_ratio = r$data_cfg$sn_ratio,
    PiLow = r$data_cfg$PiLow,
    sph = r$data_cfg$sph,
    p_signal = r$data_cfg$p_signal,
    kernel = r$method_cfg$kernel,
    criterion = r$method_cfg$criterion,
    alpha = r$method_cfg$alpha,
    fuzz_group = r$method_cfg$fuzz_group,
    fuzzifier = r$fuzzifier,
    W_eff = 1 / sum(r$W^2),
    max_W = max(r$W),
    signal_share = sum(r$W[seq_len(r$data_cfg$p_signal)]),
    ARI = r$ARI %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

# Process one seed at a time, return a data frame for that seed
process_seed <- function(seed) {
  pattern <- sprintf("^D-s%02d-", seed)
  all_files <- list.files("results",
                          pattern = "\\.rds$", full.names = TRUE)
  files <- all_files[grepl(pattern, basename(all_files))]
  if (length(files) == 0) {
    cat(sprintf("[seed %02d] no files found\n", seed))
    return(NULL)
  }
  cat(sprintf("[seed %02d] processing %d files\n", seed, length(files)))
  rows <- lapply(files, extract_one)
  do.call(rbind, rows)
}

# Parallel: one seed per core, 25 cores total
seeds <- 1:25
results <- mclapply(seeds, process_seed, mc.cores = 25, mc.preschedule = FALSE)

# Combine all seeds
weff_df <- do.call(rbind, results)
saveRDS(weff_df, "results/weff_df.rds")
cat("Saved", nrow(weff_df), "rows total across",
    length(unique(weff_df$seed)), "seeds\n")