#' Nonnegative Matrix Factorization by HALS (Frobenius loss)
#'
#' Internal single-threaded NMF used by \code{twoStageNMF()} for its first
#' stage. Minimizes \eqn{\lVert X - WH \rVert_F^2} subject to
#' \eqn{W, H \ge 0} with hierarchical alternating least squares
#' (Cichocki & Phan, 2009). Only two products touch \code{X} per iteration
#' (\code{X \%*\% t(H)} and \code{crossprod(X, W)}); the matrix \code{W \%*\% H}
#' is never formed, and the error is computed from small matrices.
#'
#' @param X A nonnegative matrix (dense, or a \code{dgCMatrix}).
#' @param rank Number of components.
#' @param max_ite Maximum number of iterations.
#' @param tol Stop when the relative change of the error is at most \code{tol}.
#' @param eps Lower bound for factor entries (keeps them strictly positive).
#'
#' @return A list with \code{W} (features x rank), \code{H} (rank x samples),
#'   \code{niter}, \code{error_out} (squared Frobenius error per iteration),
#'   \code{final_error} and \code{method = "hals"}.
#'
#' @keywords internal
#' @noRd

halsNMF <- function(X, rank, max_ite = 1000, tol = 1e-4, eps = 1e-16) {
    p <- nrow(X)
    n <- ncol(X)
    is_sparse <- inherits(X, "sparseMatrix")

    ## ||X||_F^2 and mean(X) without allocating a copy of X.
    if (is_sparse) {
        norm_X2 <- sum(X@x^2)
        mean_X <- sum(X@x) / (p * n)
    } else {
        norm_X2 <- norm(X, "F")^2
        mean_X <- sum(X) / (p * n)
    }

    ## Random start, scaled so that mean(W %*% H) matches mean(X).
    W <- matrix(runif(p * rank), nrow = p, ncol = rank)
    H <- matrix(runif(rank * n), nrow = rank, ncol = n)
    mean_WH <- sum(colSums(W) * rowSums(H)) / (p * n)
    start_scale <- sqrt(mean_X / mean_WH)
    W <- W * start_scale
    H <- H * start_scale

    error_out <- numeric(max_ite)

    for (iteration in seq_len(max_ite)) {
        ## W update: one pass over X.
        XHt <- as.matrix(X %*% t(H))          # p x rank
        HHt <- tcrossprod(H)                  # rank x rank
        for (k in seq_len(rank)) {
            W[, k] <- pmax(eps, W[, k] + (XHt[, k] - W %*% HHt[, k]) / HHt[k, k])
        }

        ## Unit-norm columns of W; H absorbs the scale (W %*% H unchanged).
        col_norm <- sqrt(colSums(W^2))
        col_norm[col_norm == 0] <- 1
        W <- sweep(W, 2L, col_norm, "/")
        H <- H * col_norm

        ## H update: one pass over X.
        WtX <- t(as.matrix(crossprod(X, W)))  # rank x n
        WtW <- crossprod(W)                   # rank x rank
        for (k in seq_len(rank)) {
            H[k, ] <- pmax(eps, H[k, ] + (WtX[k, ] - WtW[k, ] %*% H) / WtW[k, k])
        }

        ## ||X - WH||^2 = ||X||^2 - 2 <W'X, H> + <W'W, HH'>
        error_out[iteration] <- norm_X2 - 2 * sum(WtX * H) + sum(WtW * tcrossprod(H))

        if (iteration > 1L) {
            previous_error <- error_out[iteration - 1L]
            relative_change <- abs(previous_error - error_out[iteration]) /
                max(previous_error, 1e-12)
            if (relative_change <= tol) {
                break
            }
        }
    }

    error_out <- error_out[seq_len(iteration)]

    list(
        W = W,
        H = H,
        niter = iteration,
        error_out = error_out,
        final_error = error_out[iteration],
        method = "hals"
    )
}
