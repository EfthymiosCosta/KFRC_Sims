### Membership Scaling FCM (MSFCM)
.fcm_update_centers <- function(U, X, m) {
  Um <- U^m
  num <- Um %*% X
  den <- rowSums(Um)
  num / den
}

# Standard FCM membership update from distances D = (d_ij), c x n
.fcm_update_memberships <- function(D, m) {
  q <- 2 / (m - 1)
  zero_mask <- D == 0
  if (any(zero_mask)) {
    U <- matrix(0, nrow(D), ncol(D))
    for (j in which(colSums(zero_mask) > 0)) {
      zeros_j <- which(zero_mask[, j])
      U[zeros_j, j] <- 1 / length(zeros_j)
    }
    nz_cols <- which(colSums(zero_mask) == 0)
    if (length(nz_cols) > 0) {
      W <- 1 / D[, nz_cols, drop = FALSE]^q  
      U[, nz_cols] <- sweep(W, 2, colSums(W), "/")
    }
    return(U)
  }
  W <- 1 / D^q
  sweep(W, 2, colSums(W), "/")
}

# Pairwise distance matrix between cluster centres V (c x p) and data X (n x p)
.center_data_dists <- function(V, X) {
  v2 <- rowSums(V^2)
  x2 <- rowSums(X^2)
  D2 <- outer(v2, x2, "+") - 2 * V %*% t(X)
  D2[D2 < 0] <- 0 
  sqrt(D2)
}

#' Membership Scaling Fuzzy C-Means
#'
#' @param X        n x p data matrix.
#' @param c        Number of clusters.
#' @param m        Fuzziness exponent (m > 1; default 2).
#' @param eps      Convergence threshold on ||V^(t+1) - V^(t)||_F (default 1e-6).
#' @param max_iter Maximum iterations (default 200).
#' @param U_init   Optional c x n initial membership matrix. If NULL, random.
#' @param verbose  Print progress per iteration (default FALSE).
#'
#' @return A list with components:
#'   U: c x n final membership matrix
#'   V: c x p final cluster centres
#'   iter: number of iterations
#'   converged: logical
#'   J_fuzzy: fuzzy objective at final iteration (sum_i sum_j u_ij^m ||x_j - v_i||^2)
#'   J_hard: hard-clustering objective (sum_j ||x_j - v_{I*_j}||^2)
#'   filter_rate: vector of |Q_t|/n per iteration
msfcm <- function(X, c, m = 2, eps = 1e-6, max_iter = 200,
                  U_init = NULL, verbose = FALSE) {
  if (m <= 1) stop("m must be > 1")
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (c < 2 || c > n) stop("c must be between 2 and n")
  if (is.null(U_init)) {
    U <- matrix(runif(c * n), c, n)
    U <- sweep(U, 2, colSums(U), "/")
  } else {
    if (!all(dim(U_init) == c(c, n)))
      stop("U_init must be c x n")
    U <- U_init
  }
  
  V <- .fcm_update_centers(U, X, m)
  
  filter_rate <- numeric(0)
  converged <- FALSE
  
  for (t in seq_len(max_iter)) {
    
    D <- .center_data_dists(V, X) 
    U_t <- .fcm_update_memberships(D, m)
    V_tilde <- .fcm_update_centers(U_t, X, m)
    delta <- sqrt(rowSums((V_tilde - V)^2))
    delta_max <- max(delta)
    order_j  <- apply(D, 2, order)
    D_sorted <- apply(D, 2, sort)
    Istar <- order_j[1, ]
    D1 <- D_sorted[1, ]
    D2 <- D_sorted[2, ]
    Q_mask <- (D2 - delta_max) >= (D1 + delta[Istar])
    Q_idx  <- which(Q_mask)
    filter_rate <- c(filter_rate, length(Q_idx) / n)
    U_new <- U_t
    if (length(Q_idx) > 0) {
      Dc <- D_sorted[c, Q_idx]
      D1Q <- D1[Q_idx]
      IstarQ <- Istar[Q_idx]
      M <- 1 / (1 + (c - 1) * (D1Q / Dc)^(2 / (m - 1)))
      uIstar_t <- U_t[cbind(IstarQ, Q_idx)]
      beta <- ifelse(uIstar_t < 1 - 1e-12,
                     (1 - M) / (1 - uIstar_t),
                     0)
      U_new[, Q_idx] <- sweep(U_t[, Q_idx, drop = FALSE], 2, beta, "*")
      U_new[cbind(IstarQ, Q_idx)] <- M
    }
    V_new <- .fcm_update_centers(U_new, X, m)
    centre_shift <- sqrt(sum((V_new - V)^2))
    if (verbose) {
      cat(sprintf("iter %3d: filter rate = %.3f, ||dV|| = %.6g\n",
                  t, length(Q_idx) / n, centre_shift))
    }
    if (length(Q_idx) == n) {
      U <- U_new; V <- V_new
      converged <- TRUE
      break
    }
    if (centre_shift < eps) {
      U <- U_new; V <- V_new
      converged <- TRUE
      break
    }
    U <- U_new
    V <- V_new
  }
  
  D_final <- .center_data_dists(V, X)
  J_fuzzy <- sum((U^m) * D_final^2)
  hard <- max.col(t(U), ties.method = "first")
  J_hard <- sum(D_final[cbind(hard, seq_len(n))]^2)
  
  list(
    U = U,
    V = V,
    iter = t,
    converged = converged,
    J_fuzzy = J_fuzzy,
    J_hard = J_hard,
    filter_rate = filter_rate,
    hard = hard
  )
}

### Fuzzy C-Means with Varying Fuzziness Parameter (vFCM)

# FCM centre update with fuzziness parameter m
.fcm_update_centers <- function(U, X, m) {
  Um <- U^m
  num <- Um %*% X
  den <- rowSums(Um)
  num / den
}

# FCM membership update
.fcm_update_memberships <- function(D, m) {
  q <- 2 / (m - 1)
  zero_mask <- D == 0
  if (any(zero_mask)) {
    U <- matrix(0, nrow(D), ncol(D))
    for (j in which(colSums(zero_mask) > 0)) {
      zeros_j <- which(zero_mask[, j])
      U[zeros_j, j] <- 1 / length(zeros_j)
    }
    nz_cols <- which(colSums(zero_mask) == 0)
    if (length(nz_cols) > 0) {
      W <- 1 / D[, nz_cols, drop = FALSE]^q
      U[, nz_cols] <- sweep(W, 2, colSums(W), "/")
    }
    return(U)
  }
  W <- 1 / D^q
  sweep(W, 2, colSums(W), "/")
}

# Pairwise distances between centres V and data X
.center_data_dists <- function(V, X) {
  v2 <- rowSums(V^2)
  x2 <- rowSums(X^2)
  D2 <- outer(v2, x2, "+") - 2 * V %*% t(X)
  D2[D2 < 0] <- 0
  sqrt(D2)
}

#' @param X        n x p data matrix.
#' @param c        Number of clusters.
#' @param m0       Initial fuzziness exponent (default 2; authors recommend 2).
#' @param a        Multiplicative annealing factor (default 0.95).
#' @param b        Additive annealing offset (default 0.05). Together with a,
#'                 the schedule m_{t+1} = a*m_t + b has fixed point b/(1-a) = 1.
#' @param k        Update m every k iterations (default 2).
#' @param eps      Convergence threshold on ||V^(t+1) - V^(t)||_F (default 1e-6).
#' @param max_iter Maximum iterations (default 200).
#' @param U_init   Optional c x n initial membership matrix. If NULL, random.
#' @param verbose  Print progress per iteration (default FALSE).
#'
#' @return A list with components:
#'   U: c x n final membership matrix
#'   V: c x p final cluster centres
#'   m_final: final fuzziness parameter when convergence was reached
#'   iter: number of iterations
#'   converged: logical
#'   J_fuzzy: fuzzy objective at final m and final (U, V)
#'   J_hard: hard-clustering objective sum_j ||x_j - v_{I*_j}||^2
#'   m_history: vector of m values used per iteration
#'   hard: hard cluster assignments via argmax(U)
vfcm <- function(X, c, m0 = 2, a = 0.95, b = 0.05, k = 2,
                 eps = 1e-6, max_iter = 200,
                 U_init = NULL, verbose = FALSE) {
  if (m0 <= 1) stop("m0 must be > 1")
  if (a <= 0 || a >= 1) stop("a must be in (0, 1)")
  if (b <= 0) stop("b must be positive (otherwise m_t collapses below 1)")
  if (k < 1) stop("k must be >= 1")
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (c < 2 || c > n) stop("c must be between 2 and n")
  if (is.null(U_init)) {
    U <- matrix(runif(c * n), c, n)
    U <- sweep(U, 2, colSums(U), "/")
  } else {
    if (!all(dim(U_init) == c(c, n))) stop("U_init must be c x n")
    U <- U_init
  }
  m <- m0
  V <- .fcm_update_centers(U, X, m)
  m_history <- numeric(0)
  converged <- FALSE
  
  for (t in seq_len(max_iter)) {
    D <- .center_data_dists(V, X)
    U_new <- .fcm_update_memberships(D, m)
    V_new <- .fcm_update_centers(U_new, X, m)
    m_history <- c(m_history, m)
    if (t %% k == 0) {
      m <- max(a * m + b, 1.01) 
    }
    centre_shift <- sqrt(sum((V_new - V)^2))
    if (verbose) {
      cat(sprintf("iter %3d: m = %.4f, ||dV|| = %.6g\n",
                  t, m_history[t], centre_shift))
    }
    U <- U_new
    V <- V_new
    if (centre_shift < eps) {
      converged <- TRUE
      break
    }
  }
  D_final <- .center_data_dists(V, X)
  J_fuzzy <- sum((U^m_history[length(m_history)]) * D_final^2)
  hard <- max.col(t(U), ties.method = "first")
  J_hard <- sum(D_final[cbind(hard, seq_len(n))]^2)
  list(
    U = U,
    V = V,
    m_final = m_history[length(m_history)],
    iter = t,
    converged = converged,
    J_fuzzy = J_fuzzy,
    J_hard = J_hard,
    m_history = m_history,
    hard = hard
  )
}

### Fuzzy C-Means with Multiple Kernels (FCM-MK)

#' Default bandwidth grid suggested in the paper:
default_bandwidths <- function(X, multipliers = c(0.1, 0.2, 0.3, 0.4)) {
  D <- sqrt(sum((apply(X, 2, max) - apply(X, 2, min))^2))
  multipliers * D
}

# Squared distances between centres V and data X
.squared_dists <- function(V, X) {
  v2 <- rowSums(V^2)
  x2 <- rowSums(X^2)
  D2 <- outer(v2, x2, "+") - 2 * V %*% t(X)
  D2[D2 < 0] <- 0
  D2
}

# Per-kernel Gaussian similarities
.per_kernel_sims <- function(D2, sigmas) {
  S <- length(sigmas)
  arr <- array(0, c(nrow(D2), ncol(D2), S))
  for (l in seq_len(S)) {
    arr[, , l] <- exp(-D2 / (2 * sigmas[l]^2))
  }
  arr
}

# Cluster-specific normalised kernel K_tilde^(i)(x_j, v_i)
.K_tilde <- function(K_arr, W, sigmas) {
  S <- dim(K_arr)[3]
  c_n <- W / matrix(sigmas, nrow(W), S, byrow = TRUE) 
  denom <- rowSums(c_n)
  num <- matrix(0, dim(K_arr)[1], dim(K_arr)[2])
  for (l in seq_len(S)) {
    num <- num + c_n[, l] * K_arr[, , l]
  }
  num / denom
}

# Centre-update kernel K_bar^(i)(x_j, v_i)
.K_bar <- function(K_arr, W, sigmas) {
  S <- dim(K_arr)[3]
  num_coef <- W / matrix(sigmas^3, nrow(W), S, byrow = TRUE)
  den_coef <- W / matrix(sigmas, nrow(W), S, byrow = TRUE)
  denom <- rowSums(den_coef)
  num <- matrix(0, dim(K_arr)[1], dim(K_arr)[2])
  for (l in seq_len(S)) {
    num <- num + num_coef[, l] * K_arr[, , l]
  }
  num / denom
}

# Membership update: standard FCM on dist^2 = 2 - 2 K_tilde
.fcm_memberships <- function(dist2, m) {
  q <- 1 / (m - 1)
  zero_mask <- dist2 == 0
  if (any(zero_mask)) {
    U <- matrix(0, nrow(dist2), ncol(dist2))
    for (j in which(colSums(zero_mask) > 0)) {
      zeros_j <- which(zero_mask[, j])
      U[zeros_j, j] <- 1 / length(zeros_j)
    }
    nz_cols <- which(colSums(zero_mask) == 0)
    if (length(nz_cols) > 0) {
      W_inv <- 1 / dist2[, nz_cols, drop = FALSE]^q
      U[, nz_cols] <- sweep(W_inv, 2, colSums(W_inv), "/")
    }
    return(U)
  }
  W_inv <- 1 / dist2^q
  sweep(W_inv, 2, colSums(W_inv), "/")
}

# Centre update using K_bar
.fcm_centers <- function(U, X, K_bar_mat, m) {
  Um <- U^m
  WK <- Um * K_bar_mat
  num <- WK %*% X
  den <- rowSums(WK)
  num / den
}

# Compute current objective J
.objective <- function(U, K_tilde_mat, m) {
  2 * sum((U^m) * (1 - K_tilde_mat))
}

# Gradient of J w.r.t. w_il
.weight_gradient <- function(U, K_arr, W, sigmas, m) {
  S <- length(sigmas)
  c_count <- nrow(W)
  Um <- U^m
  K_tilde_mat <- .K_tilde(K_arr, W, sigmas)
  denom <- rowSums(W / matrix(sigmas, c_count, S, byrow = TRUE))
  G <- matrix(0, c_count, S)
  for (l in seq_len(S)) {
    diff_lj <- K_arr[, , l] - K_tilde_mat
    inner <- rowSums(Um * diff_lj)
    G[, l] <- -2 * inner / (sigmas[l] * denom)
  }
  G
}

# Project a vector onto the probability simplex
.project_simplex <- function(v) {
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  rho <- max(which(u - cssv / seq_along(u) > 0))
  theta <- cssv[rho] / rho
  pmax(v - theta, 0)
}

# Update W via projected gradient descent with backtracking line search.
.update_weights <- function(U, K_arr, W, sigmas, m,
                            rho_init = 0.1, rho_min = 1e-6,
                            max_inner = 50, w_eps = 1e-6) {
  S <- length(sigmas)
  c_count <- nrow(W)
  W_old <- W
  
  for (q in seq_len(max_inner)) {
    G <- .weight_gradient(U, K_arr, W, sigmas, m)
    rho <- rho_init
    accepted <- FALSE
    J_curr <- .objective(U, .K_tilde(K_arr, W, sigmas), m)
    while (rho > rho_min) {
      W_try <- W - rho * G
      # Project each row onto the simplex
      W_try <- t(apply(W_try, 1, .project_simplex))
      J_new <- .objective(U, .K_tilde(K_arr, W_try, sigmas), m)
      if (J_new < J_curr) {
        accepted <- TRUE
        break
      }
      rho <- rho / 2
    }
    if (!accepted) break
    delta <- max(abs(W_try - W))
    W <- W_try
    if (delta < w_eps) break
  }
  W
}

#' @param X         n x p data matrix.
#' @param c         Number of clusters.
#' @param sigmas    Numeric vector of kernel bandwidths. If NULL, uses
#'                  default_bandwidths(X).
#' @param m         Fuzzifier (default 2).
#' @param max_iter  Outer-loop iteration cap (default 100).
#' @param eps       Convergence tolerance on ||V^(t+1) - V^(t)||_F (default 1e-5).
#' @param rho_init  Initial gradient-descent step size for the weight update.
#' @param U_init    Optional c x n initial membership matrix. If NULL, random.
#' @param V_init    Optional c x p initial centres. If NULL, set from U_init.
#' @param verbose   Print progress per iteration (default FALSE).
#'
#' @return List with U (c x n memberships), V (c x p centres), W (c x S
#'         resolution weights), iter, converged, J (objective trace), and
#'         hard cluster assignments.
fcm_mk <- function(X, c, sigmas = NULL, m = 2,
                   max_iter = 100, eps = 1e-5,
                   rho_init = 0.1,
                   U_init = NULL, V_init = NULL,
                   verbose = FALSE) {
  if (m <= 1) stop("m must be > 1")
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (c < 2 || c > n) stop("c must be between 2 and n")
  
  if (is.null(sigmas)) sigmas <- default_bandwidths(X)
  S <- length(sigmas)
  if (is.null(U_init)) {
    U <- matrix(runif(c * n), c, n)
    U <- sweep(U, 2, colSums(U), "/")
  } else {
    if (!all(dim(U_init) == c(c, n))) stop("U_init must be c x n")
    U <- U_init
  }
  
  W <- matrix(1 / S, c, S)
  if (is.null(V_init)) {
    V <- (U^m %*% X) / rowSums(U^m)
  } else {
    if (!all(dim(V_init) == c(c, p))) stop("V_init must be c x p")
    V <- V_init
  }
  
  J_trace <- numeric(0)
  converged <- FALSE
  
  for (t in seq_len(max_iter)) {
    
    D2 <- .squared_dists(V, X)
    K_arr <- .per_kernel_sims(D2, sigmas)
    K_tilde_mat <- .K_tilde(K_arr, W, sigmas)
    K_bar_mat <- .K_bar(K_arr, W, sigmas)
    V_new <- .fcm_centers(U, X, K_bar_mat, m)
    D2_new <- .squared_dists(V_new, X)
    K_arr_new <- .per_kernel_sims(D2_new, sigmas)
    K_tilde_new <- .K_tilde(K_arr_new, W, sigmas)
    dist2 <- 2 - 2 * K_tilde_new
    dist2[dist2 < 0] <- 0
    U_new <- .fcm_memberships(dist2, m)
    W_new <- .update_weights(U_new, K_arr_new, W, sigmas, m,
                             rho_init = rho_init)
    K_tilde_final <- .K_tilde(K_arr_new, W_new, sigmas)
    J_t <- .objective(U_new, K_tilde_final, m)
    J_trace <- c(J_trace, J_t)
    centre_shift <- sqrt(sum((V_new - V)^2))
    if (verbose) {
      cat(sprintf("iter %3d: ||dV|| = %.6g, J = %.6f\n",
                  t, centre_shift, J_t))
    }
    
    U <- U_new
    V <- V_new
    W <- W_new
    if (centre_shift < eps) { converged <- TRUE; break }
  }
  
  hard <- max.col(t(U), ties.method = "first")
  list(
    U = U, V = V, W = W,
    sigmas = sigmas,
    iter = t,
    converged = converged,
    J = J_trace,
    hard = hard
  )
}