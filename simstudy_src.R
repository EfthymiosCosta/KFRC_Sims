# Kernel matrix computation

#' Compute Gaussian (RBF) kernel Gram matrix
#'
#' @param X Data matrix (n x p)
#' @param bw Vector of bandwidths (length p)
#' @param kernel Kernel type: "product" (default) or "sum"
#' @return Kernel Gram matrix K (n x n)
compute_kernel_matrix <- function(X, bw, kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (kernel == "product") {
    K <- matrix(1, n, n)
    for (j in 1:p) {
      diff_sq <- outer(X[,j], X[,j], function(a,b) (a-b)^2)
      K <- K * exp(-diff_sq / (2 * bw[j]^2))
    }
  } else {
    K <- matrix(0, n, n)
    for (j in 1:p) {
      diff_sq <- outer(X[,j], X[,j], function(a,b) (a-b)^2)
      K <- K + exp(-diff_sq / (2 * bw[j]^2))
    }
  }
  return(K)
}

# Gradient computation

#' Compute gradient of kernel matrix w.r.t. bandwidth h_j
#'
#' Product: dK/dh_j = K * diff_sq_j / h_j^3
#' Sum: dK/dh_j = K_j * diff_sq_j / h_j^3  (only jth component)
#'
#' @param X Data matrix (n x p)
#' @param K Current full kernel matrix
#' @param bw Current bandwidth vector
#' @param j Dimension index
#' @param kernel Kernel type: "product" or "sum"
#' @return Gradient matrix dK/dh_j (n x n)
gradient_K_wrt_bw <- function(X, K, bw, j, kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  diff_sq <- outer(X[,j], X[,j], function(a,b) (a-b)^2)
  if (kernel == "product") {
    dK <- K * diff_sq / bw[j]^3
  } else {
    K_j <- exp(-diff_sq / (2 * bw[j]^2))
    dK <- K_j * diff_sq / bw[j]^3
  }
  return(dK)
}

#' Compute gradient of eigenvalue i w.r.t. log-bandwidth
#'
#' @param X Data matrix
#' @param bw Current bandwidth vector
#' @param F_result Result from compute_F
#' @param i Eigenvalue index (1 = largest)
#' @param kernel Kernel type: "product" or "sum"
#' @return Gradient vector (length p)
gradient_lambda_i_wrt_log_bw <- function(X, bw, F_result, i = 1,
                                         kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  n <- nrow(X); p <- ncol(X)
  K <- compute_kernel_matrix(X, bw, kernel)
  K_bar <- F_result$K_bar
  d <- F_result$D
  v <- F_result$eigenvectors[, i]
  inv_sqrt_d <- 1 / sqrt(d)
  w <- v * inv_sqrt_d
  w_c <- w - mean(w)
  Kbar_w <- K_bar %*% w
  grad <- numeric(p)
  for (j in 1:p) {
    diff_sq <- outer(X[,j], X[,j], function(a,b) (a-b)^2)
    if (kernel == "product") {
      dK_j <- K * diff_sq / bw[j]^3
    } else {
      K_j <- exp(-diff_sq / (2 * bw[j]^2))
      dK_j <- K_j * diff_sq / bw[j]^3
    }
    dd_j <- -2 * rowMeans(dK_j) + mean(dK_j)
    dKj_wc <- dK_j %*% w_c
    dKbar_j_w <- dKj_wc - mean(dKj_wc)
    t1 <- (1/n) * as.numeric(crossprod(w, dKbar_j_w))
    t2 <- -(1/n) * as.numeric(crossprod((dd_j / d^1.5) * v, Kbar_w))
    grad[j] <- bw[j] * (t1 + t2)
  }
  return(grad)
}

#' Compute gradient of sum(log(lambda_i)) for i = 1...(c-1) w.r.t. log-bandwidth
#'
#' Gradient is sum_i (1/lambda_i) * grad(lambda_i), which automatically
#' up-weights the signal from small eigenvalues relative to large ones.
#'
#' @param X Data matrix
#' @param bw Current bandwidth vector
#' @param F_result Result from compute_F (must contain c eigenpairs)
#' @param c Number of clusters
#' @param kernel Kernel type: "product" or "sum"
#' @return Gradient vector (length p)
gradient_logsum_wrt_log_bw <- function(X, bw, F_result, c,
                                       kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  eigenvalues <- F_result$eigenvalues
  grad <- numeric(ncol(X))
  for (i in seq_len(c - 1)) {
    lam_i <- max(eigenvalues[i], 1e-10)
    grad_i <- gradient_lambda_i_wrt_log_bw(X, bw, F_result, i, kernel)
    grad <- grad + (1 / lam_i) * grad_i
  }
  return(grad)
}

# F Matrix construction
#'
#' Compute F matrix and its eigendecomposition
#'
#' @param X Data matrix (n x p)
#' @param bw Bandwidth vector
#' @param n_eigs Number of eigenpairs to compute (NULL = all)
#' @param kernel Kernel type: "product" or "sum"
#' @return List with F, K_bar, D, eigenvalues, eigenvectors
compute_F <- function(X, bw, n_eigs = NULL, kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  K <- compute_kernel_matrix(X, bw, kernel)
  n <- nrow(K)
  row_means <- rowMeans(K)
  grand_mean <- mean(row_means)
  K_bar <- K - row_means - rep(row_means, each = n) + grand_mean
  K_bar <- (K_bar + t(K_bar)) / 2
  d <- pmax(diag(K_bar), 1e-10)
  inv_sqrt_d <- 1 / sqrt(d)
  F_mat <- (1/n) * inv_sqrt_d * K_bar * rep(inv_sqrt_d, each = n)
  F_mat <- (F_mat + t(F_mat)) / 2
  if (!is.null(n_eigs) && requireNamespace("RSpectra", quietly = TRUE)) {
    eig <- RSpectra::eigs_sym(F_mat, k = n_eigs, which = "LM")
  } else {
    eig <- eigen(F_mat, symmetric = TRUE)
  }
  list(F = F_mat, K_bar = K_bar, D = d,
       eigenvalues = eig$values, eigenvectors = eig$vectors)
}

#' Compute kernel-based dissimilarity matrix
#'
#' Product kernel: R_jk = sqrt(2*(1 - K_jk)), since K_ii = 1
#' Sum kernel: R_jk = sqrt(2*(p - K_jk)), since K_ii = p
#'
#' @param X Data matrix
#' @param bw Bandwidth vector
#' @param kernel Kernel type: "product" or "sum"
#' @return Dissimilarity matrix
compute_dissimilarity <- function(X, bw, kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  K <- compute_kernel_matrix(X, bw, kernel)
  p <- ncol(X)
  diag_val <- if (kernel == "product") 1 else p
  sqrt(pmax(2 * (diag_val - K), 0))
}

#' Compute kernel-based squared dissimilarity matrix
#'
#' @param X Data matrix
#' @param bw Bandwidth vector
#' @param kernel Kernel type: "product" or "sum"
#' @return Squared dissimilarity matrix
compute_dissimilarity_sq <- function(X, bw, kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  K <- compute_kernel_matrix(X, bw, kernel)
  p <- ncol(X)
  diag_val <- if (kernel == "product") 1 else p
  2 * (diag_val - K)
}

# Projection helpers
#
#' Project log-bandwidths onto anisotropy constraint
#'
#' Ensures max(eta) - min(eta) <= log(alpha), i.e. max(bw)/min(bw) <= alpha
#'
#' @param eta Log-bandwidth vector
#' @param log_alpha Log of the anisotropy bound
#' @return Projected log-bandwidth vector
project_anisotropy <- function(eta, log_alpha) {
  spread <- max(eta) - min(eta)
  if (spread > log_alpha) {
    eta_mean <- mean(eta)
    eta <- eta_mean + (log_alpha / spread) * (eta - eta_mean)
  }
  return(eta)
}

#' Project log-bandwidths to satisfy lambda_1(F) >= tau via backtracking
#'
#' @param X Data matrix
#' @param eta Current log-bandwidth vector
#' @param eta_new Proposed new log-bandwidth vector
#' @param tau Stability threshold
#' @param step_max Maximum step factor
#' @param min_factor Minimum step factor before giving up
#' @param kernel Kernel type: "product" or "sum"
#' @return Accepted eta satisfying lambda_1 >= tau - 1e-6
project_feasibility <- function(X, eta, eta_new, tau = 0.5,
                                step_max = 1, min_factor = 1e-4,
                                kernel = c("product", "sum")) {
  kernel <- match.arg(kernel)
  step <- step_max
  while (step > min_factor) {
    eta_try <- eta + step * (eta_new - eta)
    F_try <- compute_F(X, exp(eta_try), n_eigs = 1, kernel = kernel)
    if (F_try$eigenvalues[1] >= tau - 1e-6) return(eta_try)
    step <- step / 2
  }
  warning("Could not restore lambda_1 >= tau; using original eta")
  return(eta)
}

# Two-stage bandwidth selection
#
#' Two-stage bandwidth selection for Kernel Fuzzy Relational Clustering
#'
#' Stage 1: Gradient ascent to maximise lambda_1(F), collecting all starts
#'          that achieve lambda_1 >= tau (feasible solutions).
#' Stage 2: From each feasible start, optimise a multi-cluster criterion
#'          subject to lambda_1 >= tau via backtracking projection.
#'
#' @param X Data matrix (n x p)
#' @param c Number of clusters
#' @param tau Stability threshold; require lambda_1(F) >= tau (default 0.5)
#' @param alpha Anisotropy bound: max(bw)/min(bw) <= alpha (default Inf)
#' @param kernel Kernel type: "product" (default) or "sum"
#' @param criterion Stage 2 objective: "sum", "gap", or "logsum" (default)
#' @param n_starts Number of random initialisations for Stage 1 (default 5)
#' @param max_iter Maximum gradient-ascent iterations per stage (default 200)
#' @param beta Momentum (default 0.9)
#' @param lr Learning rate (default 0.15)
#' @param tol Convergence tolerance on eta (default 1e-5)
#' @param verbose Print progress (default TRUE)
#' @return List with bw, lambda1, eigenvalues, criterion_value, eigengap,
#'         feasible, best_start_id, stage1_feasible_count, stage2_results
select_bandwidth <- function(X, c, tau = 0.5, alpha = Inf,
                             kernel = c("product", "sum"),
                             criterion = c("logsum", "sum", "gap", "minmax"),
                             n_starts = 5, max_iter = 200, beta = 0.9,
                             lr = 0.15, tol = 1e-5, verbose = TRUE) {
  kernel <- match.arg(kernel, c("product", "sum"))
  criterion <- match.arg(criterion, c("logsum", "sum", "gap", "minmax"))
  X <- as.matrix(X)
  p <- ncol(X)
  log_alpha <- log(alpha)
  sds <- apply(X, 2, sd)
  if (verbose) {
    cat("Stage 1: Maximising lambda_1 with", n_starts, "random starts\n")
    cat(sprintf("  tau = %.2f,  kernel = '%s'", tau, kernel))
    if (is.finite(alpha)) cat(sprintf(",  alpha = %.2f", alpha))
    cat("\n")
  }
  feasible_starts <- list()
  for (start in seq_len(n_starts)) {
    eta <- if (start == 1) log(sds) else log(sds) + rnorm(p, 0, 0.5)
    if (is.finite(alpha)) eta <- project_anisotropy(eta, log_alpha)
    v <- numeric(length(eta))
    for (iter in seq_len(max_iter)) {
      bw <- exp(eta)
      F_res <- compute_F(X, bw, n_eigs = 1, kernel = kernel)
      grad <- gradient_lambda_i_wrt_log_bw(X, bw, F_res, i = 1, kernel = kernel)
      v <- beta * v + lr * grad
      eta_new <- eta + v
      if (is.finite(alpha)) eta_new <- project_anisotropy(eta_new, log_alpha)
      converged <- (max(abs(eta_new - eta)) < tol)
      eta <- eta_new
      if (converged) break
    }
    bw <- exp(eta)
    F_res <- compute_F(X, bw, n_eigs = c, kernel = kernel)
    lam1 <- F_res$eigenvalues[1]
    if (lam1 >= tau) {
      feasible_starts <- append(feasible_starts, list(
        list(eta = eta, F_result = F_res, lambda1 = lam1, start_id = start)
      ))
      if (verbose) cat(sprintf("  Start %d: feasible (lambda_1 = %.4f)\n", start, lam1))
    } else {
      if (verbose) cat(sprintf("  Start %d: not feasible (lambda_1 = %.4f)\n", start, lam1))
    }
  }
  if (length(feasible_starts) == 0) {
    warning("No feasible solution found in Stage 1.")
    return(list(bw = NULL, lambda1 = NULL, feasible = FALSE,
                eigenvalues = NULL, stage1_feasible_count = 0L))
  }
  if (verbose) cat(sprintf("Stage 1 complete: %d feasible solution(s)\n", length(feasible_starts)))
  if (verbose) {
    cat(sprintf("\nStage 2: Optimising criterion = '%s' (c = %d) from each feasible start\n",
                criterion, c))
  }
  stage2_results <- vector("list", length(feasible_starts))
  for (i in seq_along(feasible_starts)) {
    s_info <- feasible_starts[[i]]
    eta <- s_info$eta
    F_res <- s_info$F_result
    if (verbose) cat(sprintf("\n  Start %d (lambda_1 = %.4f)...\n",
                             s_info$start_id, s_info$lambda1))
    for (iter in seq_len(max_iter)) {
      bw <- exp(eta)
      grad_crit <- if (criterion == "sum") {
        Reduce("+", lapply(seq_len(c - 1), function(k)
          gradient_lambda_i_wrt_log_bw(X, bw, F_res, k, kernel = kernel)))
      } else if (criterion == "logsum") {
        gradient_logsum_wrt_log_bw(X, bw, F_res, c, kernel = kernel)
      } else if (criterion == "gap") {
        gradient_lambda_i_wrt_log_bw(X, bw, F_res, c - 1, kernel = kernel) -
          gradient_lambda_i_wrt_log_bw(X, bw, F_res, c,     kernel = kernel)
      } else {
        gradient_lambda_i_wrt_log_bw(X, bw, F_res, c - 1, kernel = kernel)
      }
      eta_new <- eta + (lr / 2) * grad_crit
      if (is.finite(alpha)) eta_new <- project_anisotropy(eta_new, log_alpha)
      eta_new <- project_feasibility(X, eta, eta_new, tau = tau, kernel = kernel)
      converged <- max(abs(eta_new - eta)) < tol
      eta <- eta_new
      F_res <- compute_F(X, exp(eta), n_eigs = c, kernel = kernel)
      if (converged) break
    }
    eigs_final <- F_res$eigenvalues
    crit_val <- if (criterion == "sum") {
      sum(eigs_final[seq_len(c - 1)])
    } else if (criterion == "logsum") {
      sum(log(pmax(eigs_final[seq_len(c - 1)], 1e-10)))
    } else if (criterion == "gap") {
      eigs_final[c - 1] - eigs_final[c]
    } else {
      eigs_final[c - 1]
    }
    stage2_results[[i]] <- list(
      start_id = s_info$start_id,
      bw = exp(eta),
      eta = eta,
      lambda1 = eigs_final[1],
      eigenvalues = eigs_final,
      criterion_value = crit_val,
      eigengap = if (c >= 2) eigs_final[c - 1] - eigs_final[c] else NA,
      iterations = iter,
      feasible = eigs_final[1] >= tau
    )
    if (verbose) {
      crit_label <- if (criterion == "sum") sprintf("sum(1:%d)", c - 1) else
        if (criterion == "logsum") sprintf("logsum(1:%d)", c - 1) else
          if (criterion == "gap") "eigengap" else
            sprintf("lambda_%d", c - 1)
      cat(sprintf("lambda_1 = %.4f,  %s = %.4f\n",
                  eigs_final[1], crit_label, crit_val))
    }
  }
  crit_values <- sapply(stage2_results, `[[`, "criterion_value")
  best_idx <- which.max(crit_values)
  best <- stage2_results[[best_idx]]
  if (verbose) {
    crit_label <- if (criterion == "sum") sprintf("sum(1:%d)", c - 1) else
      if (criterion == "logsum") sprintf("logsum(1:%d)", c - 1) else
        if (criterion == "gap") "eigengap" else
          sprintf("lambda_%d", c - 1)
    cat(sprintf("\nBest result: start %d | %s = %.4f | lambda_1 = %.4f\n",
                best$start_id, crit_label, best$criterion_value, best$lambda1))
    cat(sprintf("Bandwidths: %s\n", paste(round(best$bw, 3), collapse = ", ")))
    if (is.finite(alpha))
      cat(sprintf("Anisotropy ratio: %.3f (bound: %.1f)\n",
                  max(best$bw) / min(best$bw), alpha))
  }
  start_ids <- sapply(feasible_starts, `[[`, "start_id")
  stage1_bw <- setNames(lapply(feasible_starts, function(s) exp(s$eta)), start_ids)
  stage1_lambda1 <- setNames(sapply(feasible_starts, `[[`, "lambda1"), start_ids)
  list(
    bw = best$bw,
    lambda1 = best$lambda1,
    eigenvalues = best$eigenvalues,
    criterion_value = best$criterion_value,
    eigengap = best$eigengap,
    feasible = best$feasible,
    best_start_id = best$start_id,
    stage1_feasible_count = length(feasible_starts),
    stage1_bw = stage1_bw,
    stage1_lambda1 = stage1_lambda1,
    stage2_results = stage2_results
  )
}

#' Fuzzy Relational Clustering (FRC)
#'
#' @param Dmat Dissimilarity matrix (n x n)
#' @param k Number of clusters
#' @param m Fuzzifier parameter (m > 1)
#' @param maxit Maximum number of iterations
#' @param conv Convergence tolerance
#' @param U_init Optional initial membership matrix (n x k)
#'
#' @return List with U (membership matrix), iter, converged, objective
FRC <- function(Dmat, k, m = 2, maxit = 1000, conv = 1e-6, U_init = NULL) {
  n <- nrow(Dmat)
  
  if (is.null(U_init)) {
    U <- matrix(runif(n * k), n, k)
    U <- U / rowSums(U)
  } else {
    U <- U_init
  }
  
  for (iter in 1:maxit) {
    U_old <- U
    Um <- U^m
    S <- colSums(Um)
    
    Q <- numeric(k)
    for (i in 1:k) {
      Q[i] <- sum(Um[, i] * (Dmat %*% Um[, i]))
    }
    
    a <- matrix(0, nrow = n, ncol = k)
    for (i in 1:k) {
      T_i <- as.numeric(Dmat %*% Um[, i])
      a[, i] <- m * T_i / S[i] - m * Q[i] / (2 * S[i]^2)
    }
    
    U <- matrix(0, nrow = n, ncol = k)
    for (j in 1:n) {
      neg_idx <- which(a[j, ] <= 0)
      if (length(neg_idx) > 0) {
        U[j, neg_idx] <- 1 / length(neg_idx)
      } else {
        inv_a_pow <- (1 / a[j, ])^(1 / (m - 1))
        U[j, ] <- inv_a_pow / sum(inv_a_pow)
      }
    }
    
    if (max(abs(U - U_old)) < conv) break
  }
  
  Um <- U^m
  S <- colSums(Um)
  obj <- 0
  for (i in 1:k) {
    obj <- obj + sum(Um[, i] * (Dmat %*% Um[, i])) / (2 * S[i])
  }
  
  list(
    U = U,
    iter = iter,
    converged = iter < maxit,
    objective = obj
  )
}

#' Exponential fuzzifier function
#' t(u) = (exp(m*u) - 1) / (exp(m) - 1)
#'
#' @param u Membership value(s) in [0, 1]
#' @param m Fuzzifier parameter (m > 0)
t_exp <- function(u, m) {
  (exp(m * u) - 1) / (exp(m) - 1)
}

#' Fuzzy Relational Clustering with Exponential Fuzzifier (FRC_exp)
#'
#' @param Dmat Dissimilarity matrix (n x n)
#' @param k Number of clusters
#' @param m Fuzzifier parameter (m > 0)
#' @param maxit Maximum number of iterations
#' @param conv Convergence tolerance
#' @param U_init Optional initial membership matrix (n x k)
#'
#' @return List with U (membership matrix), iter, converged, objective
FRC_exp <- function(Dmat, k, m = 2, maxit = 1000, conv = 1e-6, U_init = NULL) {
  stopifnot(m > 0)
  Dmat <- as.matrix(Dmat)
  n <- nrow(Dmat)
  k <- as.integer(k)
  
  if (is.null(U_init)) {
    U <- matrix(runif(n * k), n, k)
    U <- U / rowSums(U)
  } else {
    U <- U_init
  }
  
  for (iter in seq_len(maxit)) {
    U_old <- U
    tU <- t_exp(U, m)
    S <- colSums(tU)
    RtU <- Dmat %*% tU
    Q <- colSums(tU * RtU)
    
    A <- sweep(2 * RtU, 2, S, "/") - matrix(Q / S^2, nrow = n, ncol = k, byrow = TRUE)
    A <- pmax(A, 1e-10)
    
    log_A <- log(A)        
    row_mean_lA <- rowMeans(log_A) 
    U <- 1/k + (1/m) * (row_mean_lA - log_A) 
    
    U <- pmax(U, 0)
    rs <- rowSums(U)
    zero_rows <- rs == 0
    U[!zero_rows, ] <- U[!zero_rows, ] / rs[!zero_rows]
    U[ zero_rows, ] <- 1/k
    
    if (max(abs(U - U_old)) < conv) break
  }
  
  tU <- t_exp(U, m)
  S <- colSums(tU)
  obj <- sum(colSums(tU * (Dmat %*% tU)) / (2 * S))
  
  list(U = U, iter = iter, converged = iter < maxit, objective = obj)
}

#' Complementary root fuzzifier function
#' t(u) = 1 - (1-u)^{1/m}
#' #'
#' @param u Membership value(s) in [0, 1]
#' @param m Fuzzifier parameter (m > 1)
t_rat2 <- function(u, m) 1 - (1 - u)^(1/m)

#' Fuzzy Relational Clustering with Complementary Root Fuzzifier (FRC_rd)
#'
#' @param Dmat Dissimilarity matrix (n x n)
#' @param k Number of clusters
#' @param m Fuzzifier parameter (m > 0)
#' @param maxit Maximum number of iterations
#' @param conv Convergence tolerance
#' @param U_init Optional initial membership matrix (n x k)
#'
#' @return List with U (membership matrix), iter, converged, objective
FRC_rc <- function(Dmat, k, m = 2, maxit = 1000, conv = 1e-6,
                   U_init = NULL) {
  stopifnot(m > 1)
  Dmat <- as.matrix(Dmat)
  n <- nrow(Dmat)
  k <- as.integer(k)
  e <- m / (m - 1)
  
  if (is.null(U_init)) {
    U <- matrix(runif(n * k), n, k)
    U <- U / rowSums(U)
  } else {
    U <- U_init
  }
  
  j_seq <- matrix(seq_len(k), nrow = n, ncol = k, byrow = TRUE)
  
  for (iter in seq_len(maxit)) {
    U_old <- U
    tU <- t_rat2(U, m)
    S <- colSums(tU)
    RtU <- Dmat %*% tU
    Q <- colSums(tU * RtU)
    
    A <- sweep(2 * RtU, 2, S, "/") -
      matrix(Q / S^2, nrow = n, ncol = k, byrow = TRUE)
    A <- pmax(A, 1e-10)
    
    Ae <- A^e
    sorted_idx <- t(apply(Ae, 1, order))
    Ae_sorted <- t(apply(Ae, 1, sort))
    cs <- t(apply(Ae_sorted, 1, cumsum))
    
    cond <- (j_seq - 1) * Ae_sorted < cs
    J_active <- apply(cond, 1, function(x) max(which(x)))
    
    active_mask <- j_seq <= matrix(J_active, nrow = n, ncol = k)
    cs_at_J <- cs[cbind(seq_len(n), J_active)]
    
    U_sorted <- active_mask *
      (1 - (matrix(J_active, nrow = n, ncol = k) - 1) *
         Ae_sorted / matrix(cs_at_J, nrow = n, ncol = k))
    
    U <- matrix(0, nrow = n, ncol = k)
    U[cbind(rep(seq_len(n), each = k),
            as.vector(t(sorted_idx)))] <- as.vector(t(U_sorted))
    
    if (max(abs(U - U_old)) < conv) break
  }
  
  tU <- t_rat2(U, m)
  S <- colSums(tU)
  obj <- sum(colSums(tU * (Dmat %*% tU)) / (2 * S))
  list(U = U, iter = iter, converged = iter < maxit, objective = obj)
}

#' Quadratic fuzzifier function
#'
#' @param u Membership value(s) in [0, 1]
#' @param m Fuzzifier parameter (m in [0, 1])
t_quad <- function(u, m) {
  m * u^2 + (1 - m) * u
}

#' Fuzzy Relational Clustering with Quadratic Fuzzifier (FRC_quad)
#'
#' @param Dmat Dissimilarity matrix (n x n)
#' @param k Number of clusters
#' @param m Fuzzifier parameter (m in [0, 1])
#' @param maxit Maximum number of iterations
#' @param conv Convergence tolerance
#' @param U_init Optional initial membership matrix (n x k)
#'
#' @return List with U (membership matrix), iter, converged, objective
FRC_quad <- function(Dmat, k, m = 0.5, maxit = 1000, conv = 1e-6,
                     U_init = NULL) {
  stopifnot(m > 0, m <= 1)
  Dmat <- as.matrix(Dmat)
  n <- nrow(Dmat)
  k <- as.integer(k)
  
  if (is.null(U_init)) {
    U <- matrix(runif(n * k), n, k)
    U <- U / rowSums(U)
  } else {
    U <- U_init
  }
  numerator <- k + m * (2 - k)
  for (iter in seq_len(maxit)) {
    U_old <- U
    tU <- t_quad(U, m)
    S <- colSums(tU)
    RtU <- Dmat %*% tU
    Q <- colSums(tU * RtU)
    
    A <- sweep(2 * RtU, 2, S, "/") - matrix(Q / S^2, nrow = n, ncol = k, byrow = TRUE)
    A <- pmax(A, 1e-10)
    
    inv_A <- 1 / A
    sum_inv_A <- rowSums(inv_A)
    U <- (1 / (2 * m)) * (numerator * inv_A / sum_inv_A - (1 - m))
    U <- pmax(U, 0)
    rs <- rowSums(U)
    zero_rows <- rs == 0
    U[!zero_rows, ] <- U[!zero_rows, ] / rs[!zero_rows]
    U[ zero_rows, ] <- 1/k
    if (max(abs(U - U_old)) < conv) break
  }
  tU <- t_quad(U, m)
  S <- colSums(tU)
  obj <- sum(colSums(tU * (Dmat %*% tU)) / (2 * S))
  list(U = U, iter = iter, converged = iter < maxit, objective = obj)
}

#' Fuzzy ARI (FARI_fuzzy)
#'
#' @param a A matrix of cluster memberships
#' @param b A matrix of cluster memberships
fari <- function(a, b){
  n <- nrow(a)
  A <- a %*% t(a)
  B <- b %*% t(b)
  j <- matrix(1, n, n)
  Na <- (sum(A*j)/sum(A*A))*A
  Nb <- (sum(B*j)/sum(B*B))*B
  ri <- (sum(Na*Nb) + sum((j-Na)*(j-Nb)) - n)/(2*choose(n,2))
  M <- j/n
  R <- diag(n) - M
  Eri <- ((2*sum(A*j)*sum(B*j)/(sum(A*A)*sum(B*B))) * (sum(M*A)*sum(M*B)+((1/(n-1))*sum(R*A)*sum(R*B))) - (sum(A*j)^2)/sum(A*A) - (sum(B*j)^2)/sum(B*B) + n^2 - n)/(2*choose(n,2))
  ari <- (ri - Eri)/(1 - Eri)
  ari
}

#' True posterior values for GMM
#'
#' @param X_raw_signal Raw data set X
#' @param Pi Mixture weights
#' @param Mu Mean vector
#' @param S Covariance matrix
compute_true_posteriors <- function(X_raw_signal, Pi, Mu, S) {
  n <- nrow(X_raw_signal)
  K <- length(Pi)
  log_dens <- matrix(0, n, K)
  for (k in seq_len(K))
    log_dens[, k] <- mvtnorm::dmvnorm(X_raw_signal, mean = Mu[k, ],
                                      sigma = S[, , k], log = TRUE) +
    log(Pi[k])
  max_log <- apply(log_dens, 1, max)
  log_norm <- max_log + log(rowSums(exp(log_dens - max_log)))
  exp(log_dens - log_norm)
}