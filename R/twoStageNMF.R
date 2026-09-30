#' Two-staged Nonnegative Matrix Factorization
#'
#' Two-stage decomposition of multiple matrices with Nonnegative Matrix Factorization (NMF)
#'
#' @param dataset A list of datasets to be analyzed
#' @param group A list of grouping of the datasets, indicating the relationship
#'   between the datasets.
#' @param comp_num A vector indicating the dimension of each component
#' @param comp_num_first_stage An integer indicating the number of components to
#'   be used in the first stage NMF. Must be greater or equal to the sum of
#'   `comp_num` elements associated with groups of multiple datasets.
#' @param nrun_first_stage Number of NMF runs to perform in the first stage.
#'   Multiple runs achieve greater stability and decrease chances of bad local
#'   minima, but significantly increase computational time.
#' @param proj_dataset An optional nonnegative matrix containing an additional
#'   dataset to project onto the learned second-stage components. Features must
#'   be in rows and samples in columns. If row names are available, features are
#'   matched to the learned components by name. If row names are absent, the
#'   rows are assumed to have the same order as the input datasets.
#' @param proj_group An optional logical vector of length \code{length(group)}
#'   indicating which second-stage component groups should be used to project
#'   \code{proj_dataset}. Required when \code{proj_dataset} is provided. Groups
#'   marked \code{FALSE} receive structural-zero score matrices in
#'   \code{proj_score_list}.
#' @param max_ite The maximum number of iterations for the twoStageNMF
#'   algorithms to run on its second stage for each group, default value is set
#'   to 100.
#' @param max_err The maximum error of loss between two iterations, or the
#'   program will terminate and return, default value is set to be 0.0001.
#' @param enable_normalization An argument to decide whether to use
#'   normalization or not,  default is TRUE.
#' @param column_sum_normalization An argument to decide whether to use column
#'   sum normalization or not, default it FALSE.
#' @param use_nnls Logical. If TRUE, fixed-basis projections are estimated using
#'   nonnegative least squares via `nnls::nnls()`.If FALSE, multiplicative
#'   updates are used. This does not affect sequential steps in which a new
#'   component matrix is learned. default is FALSE.
#'
#' @return A named list containing the fitted components, scores, convergence
#'   information, and intermediate results:
#'
#'   \describe{
#'
#'     \item{\code{linked_component_list}}{
#'       A list of length \code{length(group)} containing the second-stage
#'       nonnegative component matrices. Element \code{i} corresponds to
#'       \code{group[[i]]} and is a matrix with genes or features in rows and
#'       \code{comp_num[i]} components in columns. These matrices correspond
#'       to the group-specific loading matrices \eqn{G_g}. Row names are
#'       inherited from the input features, and component and list names are
#'       assigned using the SJD naming conventions.
#'     }
#'
#'     \item{\code{score_list}}{
#'       A nested list containing the scores of every original dataset on every
#'       second-stage component matrix. The outer list has one element per
#'       dataset, and the inner list has one element per group. Thus,
#'       \code{score_list[[i]][[g]]} is a matrix with
#'       \code{comp_num[g]} rows and \code{ncol(dataset[[i]])} columns.
#'
#'       If dataset \code{i} belongs to \code{group[[g]]}, the matrix contains
#'       its estimated nonnegative scores on the components associated with
#'       group \code{g}. If dataset \code{i} does not belong to that group, the
#'       corresponding matrix contains structural zeros.
#'
#'       When a singleton group is supplied for dataset \code{i}, its scores
#'       are taken from the sequential singleton fit. When no singleton group
#'       is supplied, the relevant second-stage component matrices are held
#'       fixed and projected back onto the original dataset. This projection
#'       uses multiplicative updates when \code{use_nnls = FALSE} and
#'       nonnegative least squares when \code{use_nnls = TRUE}. Dataset,
#'       group, component, and sample names follow the SJD naming conventions.
#'     }
#'
#'     \item{\code{proj_score_list}}{
#'       If \code{proj_dataset} is provided, a list with one element per
#'       second-stage group containing the scores of the projected samples.
#'       Element \code{g} has \code{comp_num[g]} rows and
#'       \code{ncol(proj_dataset)} columns. Groups selected by
#'       \code{proj_group} contain estimated scores; groups not selected contain
#'       structural-zero matrices. Scores are estimated with nonnegative least
#'       squares when \code{use_nnls = TRUE} and multiplicative updates
#'       otherwise. If \code{proj_dataset} is \code{NULL}, this element is
#'       \code{NULL}.
#'     }
#'
#'     \item{\code{proj_error_out}}{
#'       Projection error information for \code{proj_dataset}. When NNLS is
#'       used, this is the final residual sum of squares. With multiplicative
#'       updates, it is the reconstruction-error trajectory. It is \code{NULL}
#'       when no projection dataset is supplied.
#'     }
#'
#'     \item{\code{proj_method}}{
#'       The method used to estimate \code{proj_score_list}: \code{"nnls"} or
#'       \code{"multiplicative_update"}. It is \code{NULL} when no projection
#'       dataset is supplied.
#'     }
#'
#'     \item{\code{error_out}}{
#'       A list of length \code{length(group)}. Element \code{g} contains the
#'       reconstruction-error trajectory for the sequential second-stage fit
#'       of \code{group[[g]]}. Each entry records
#'       \eqn{\lVert X_g - G S \rVert_F^2} after an iteration of the
#'       multiplicative-update algorithm. The length of each vector may be
#'       smaller than \code{max_iter} when the relative change in error reaches
#'       the specified convergence tolerance.
#'     }
#'
#'     \item{\code{first_stage}}{
#'       A named list containing results from the separate first-stage NMF fits:
#'
#'       \describe{
#'
#'         \item{\code{component_list}}{
#'           A list with one element per input dataset. Element \code{i}
#'           contains the first-stage gene or feature loading matrix
#'           \eqn{G_i}, with features in rows and
#'           \code{comp_num_first_stage[i]} components in columns.
#'         }
#'
#'         \item{\code{score_list}}{
#'           A list with one element per input dataset. Element \code{i}
#'           contains the corresponding first-stage score matrix \eqn{S_i},
#'           with \code{comp_num_first_stage[i]} components in rows and the
#'           samples from dataset \code{i} in columns.
#'         }
#'
#'         \item{\code{fit_list}}{
#'           A list containing the complete fitted first-stage NMF objects
#'           returned by \code{\link[NMF]{nmf}}, one for each input dataset.
#'           These objects retain additional diagnostics and algorithm-specific
#'           information produced by the NMF package.
#'         }
#'       }
#'     }
#'
#'     \item{\code{projection_error_out}}{
#'       A list with one element per input dataset. If a dataset has no
#'       singleton group, its entry contains diagnostic information from
#'       projecting the fixed second-stage components onto that dataset.
#'
#'       For multiplicative-update projection, the entry is a numeric vector
#'       containing the reconstruction error at each iteration. For NNLS
#'       projection, the entry is a single numeric value giving the final
#'       residual sum of squares,
#'       \eqn{\lVert D_i - G_i^{*} S_i^{*} \rVert_F^2}. The entry is
#'       \code{NULL} when the dataset has a singleton group because its scores
#'       were estimated during the sequential second-stage fit rather than
#'       through a separate fixed-basis projection.
#'     }
#'
#'     \item{\code{projection_method}}{
#'       A character list with one element per input dataset indicating how its
#'       final scores were obtained. Possible values are
#'       \code{"singleton_fit"}, \code{"multiplicative_update"}, and
#'       \code{"nnls"}.
#'     }
#'
#'     \item{\code{mask_list}}{
#'       A list of structural masks used during the sequential second-stage
#'       fits. Element \code{g} has the same dimensions as the coefficient
#'       matrix used while fitting \code{group[[g]]}. Entries equal to one
#'       indicate permitted component--dataset relationships, while entries
#'       equal to zero identify coefficients constrained to remain exactly
#'       zero because the corresponding dataset does not belong to the
#'       component's group.
#'     }
#'
#'     \item{\code{target_info}}{
#'       A list with one element per second-stage group containing bookkeeping
#'       information for the corresponding sequential fit. Each element
#'       contains:
#'
#'       \describe{
#'
#'         \item{\code{group}}{
#'           The indices of the datasets included in the current target.
#'         }
#'
#'         \item{\code{column_dataset}}{
#'           An integer vector identifying the source dataset associated with
#'           each column of the target matrix used at that step.
#'         }
#'
#'         \item{\code{previous_groups}}{
#'           The indices of previously learned groups whose component matrices
#'           were included as fixed components in the current fit.
#'         }
#'       }
#'     }
#'   }
#'
#'   The first four elements,
#'   \code{linked_component_list}, \code{score_list},
#'   \code{proj_score_list}, and \code{error_out}, follow the primary output
#'   naming convention used by \code{\link{jointNMF}}. The remaining elements
#'   contain diagnostics and intermediate results specific to the two-stage
#'   procedure.
#'
#' @keywords two-staged, NMF
#'
#' @examples
#' dataset = list(matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50),
#' matrix(runif(5000, 1, 2), nrow = 100, ncol = 50))
#' group = list(c(1,2,3,4), c(1,2), c(3,4), c(1,3), c(2,4), c(1), c(2), c(3), c(4))
#' comp_num = c(2,2,2,2,2,2,2,2,2)
#' comp_num_first_stage = 30
#' res_twoStageNMF = twoStageNMF(
#' dataset,
#' group,
#' comp_num,
#' comp_num_first_stage)
#' @export

twoStageNMF <- function(dataset, group, comp_num, comp_num_first_stage, nrun_first_stage = 1, proj_dataset = NULL, proj_group = NULL, max_ite = 100, max_err = 0.0001, enable_normalization = TRUE, column_sum_normalization = FALSE, use_nnls = FALSE) {

    ## Input handling
    min_comp_first_stage = sum(comp_num[lengths(group) > 1])
    if (comp_num_first_stage < min_comp_first_stage) {
        stop("comp_num_first_stage set too low for given comp_num. \n",
        "The number of first stage dimensions must at least be equal to the sum of comp_num entries corresponding to multi-dataset groups.")
    }

    if (length(group[lengths(group) > 1]) == 0){
        stop("No group contains multiple datasets. Two-stage NMF is not appropriate in this case. Run sepNMF")
    }

    if (is.unsorted(-vapply(group, length, integer(1)))) {
        stop("`group` must be ordered from largest to smallest groups.")
    }


    ## Obtain names for dataset, gene and samples
    dataset_name = datasetNameExtractor(dataset)
    gene_name = geneNameExtractor(dataset)
    sample_name = sampleNameExtractor(dataset)
    group_name = groupNameExtractor(group)

    ## Preprocess Dataset
    dataset = frameToMatrix(dataset)

    dataset = normalizeData(dataset, enable_normalization, column_sum_normalization, nonnegative_normalization = TRUE)
    N = length(dataset)

    ## First stage NMF: separately apply NMF to each dataset and save G matrices
    first_stage_G   <- vector("list", N)
    first_stage_S   <- vector("list", N)
    first_stage_fit   <- vector("list", N)

    for (i in 1:N) {
        nmf_temp = nmf(dataset[[i]], comp_num_first_stage, nrun = nrun_first_stage)
        first_stage_fit[[i]] = nmf_temp
        first_stage_G[[i]] = nmf_temp@fit@W
        first_stage_S[[i]] = nmf_temp@fit@H
    }

    ## Second stage NMF: apply NMF to combined G matrices from first step
    K = length(group)
    p = nrow(dataset[[1]])
    first_stage_G = normalizeData(first_stage_G, enable_normalization = TRUE, column_sum_normalization = TRUE, nonnegative_normalization = TRUE)

    second_stage_G     <- vector("list", K)
    second_stage_S     <- vector("list", K)
    second_stage_mask  <- vector("list", K)
    second_stage_error <- vector("list", K)
    second_stage_target <- vector("list", K)

    for (k in 1:K) {
        current_group = group[[k]]

        ## For singleton groups, use the original dataset.
        ## Otherwise, use the concatenated first-stage bases.
        if (length(current_group) == 1L) {
            block_list <- list(dataset[[current_group]])
        } else {
            block_list <- first_stage_G[current_group]
        }

        X <- do.call(cbind, block_list)

        ## Dataset membership of every column of X.
        column_dataset <- rep(
            current_group,
            times = vapply(block_list, ncol, integer(1))
        )

        ## Previously learned groups that could explain at least one
        ## block in the current target.
        previous_groups <- if (k == 1L) {
            integer(0)
        } else {
            which(vapply(
                group[seq_len(k - 1L)],
                function(old_group) {
                    length(intersect(old_group, current_group)) > 0L
                },
                logical(1)
            ))
        }

        old_W <- if (length(previous_groups) == 0L) {
            matrix(0, nrow = p, ncol = 0L)
        } else {
            do.call(cbind, second_stage_G[previous_groups])
        }

        old_rank <- ncol(old_W)
        new_rank <- comp_num[k]

        ## Initialize the new basis.
        new_W <- matrix(
            runif(p * new_rank, min = min(X), max = max(X)),
            nrow = p,
            ncol = new_rank
        )

        W <- cbind(old_W, new_W)

        ## Build the structural-zero mask.
        ##
        ## mask[a,j] = 1 if component a may contribute to target
        ## column j, and 0 otherwise.
        mask <- matrix(
            0,
            nrow = old_rank + new_rank,
            ncol = ncol(X)
        )

        row_start <- 1L

        for (old_group in previous_groups) {
            rows <- row_start:(row_start + comp_num[old_group] - 1L)

            mask[rows, ] <- matrix(
                as.numeric(column_dataset %in% group[[old_group]]),
                nrow = length(rows),
                ncol = ncol(X),
                byrow = TRUE
            )

            row_start <- row_start + comp_num[old_group]
        }

        ## The newly introduced group applies to all columns in X,
        ## since X contains only datasets belonging to the current group.
        new_rows <- (old_rank + 1L):(old_rank + new_rank)
        mask[new_rows, ] <- 1

        H <- matrix(
            runif(nrow(mask) * ncol(X), min = min(X), max = max(X)),
            nrow = nrow(mask),
            ncol = ncol(X)
        )

        H <- H * mask

        ## Iteratively estimate the NMF with Euclidean distance
        error_history <- numeric(max_ite)

        for (iteration in seq_len(max_ite)) {
            ## Update all permitted coefficients.
            numerator_H <- crossprod(W, X)
            denominator_H <- crossprod(W, W %*% H)

            H <- H * numerator_H / denominator_H

            ## Reimpose exact structural zeros.
            H <- H * mask

            ## Update only the newly introduced basis.
            H_new <- H[new_rows, , drop = FALSE]

            numerator_W <- X %*% t(H_new)
            denominator_W <- W %*% H %*% t(H_new)

            new_W <- new_W * numerator_W / denominator_W
            new_W[!is.finite(new_W)] <- 0

            W[, new_rows] <- new_W

            error_history[iteration] <- sum((X - W %*% H)^2)

            if (iteration > 1L) {
                relative_change <- abs(
                    error_history[iteration] -
                        error_history[iteration - 1L]
                ) / max(error_history[iteration - 1L], 1e-12)

                if (relative_change <= max_err) {
                    error_history <- error_history[seq_len(iteration)]
                    break
                }
            }
        }


        second_stage_G[[k]] <- new_W
        second_stage_S[[k]] <- H
        second_stage_mask[[k]] <- mask
        second_stage_error[[k]] <- error_history

    }


    ## Reorganize output using jointNMF conventions

    linked_component_list <- second_stage_G

    ## Assign the same names used by jointNMF:
    ##   - list elements named by group
    ##   - component columns named by group
    ##   - feature/gene row names restored
    linked_component_list <- compNameAssign(
        linked_component_list,
        group_name
    )

    linked_component_list <- geneNameAssign(
        linked_component_list,
        gene_name
    )


    ## Construct score_list[[dataset]][[group]]

    estimateScoresFixedG <- function(
        X,
        G,
        use_nnls = FALSE,
        max_ite = 1000,
        max_err = 1e-5
    ) {
        X <- as.matrix(X)
        G <- as.matrix(G)

        if (nrow(X) != nrow(G)) {
            stop("`X` and `G` must have the same number of rows.")
        }

        if (ncol(G) == 0L) {
            stop("Cannot estimate scores from an empty component matrix.")
        }

        if (any(!is.finite(X)) || any(!is.finite(G))) {
            stop("`X` and `G` must contain only finite values.")
        }

        if (any(X < 0) || any(G < 0)) {
            stop("NNMF projection requires nonnegative `X` and `G`.")
        }

        if (use_nnls) {
            ## ----------------------------------------------------
            ## Exact nonnegative least-squares projection
            ##
            ## One independent NNLS problem is solved per sample:
            ##     min || X[, j] - G %*% S[, j] ||^2
            ##     subject to S[, j] >= 0
            ## ----------------------------------------------------

            S <- matrix(
                0,
                nrow = ncol(G),
                ncol = ncol(X)
            )

            for (sample_index in seq_len(ncol(X))) {
                fit <- nnls::nnls(
                    A = G,
                    b = X[, sample_index]
                )

                S[, sample_index] <- stats::coef(fit)
            }

            rownames(S) <- colnames(G)
            colnames(S) <- colnames(X)

            reconstruction_error <- sum(
                (X - G %*% S)^2
            )

            return(list(
                score = S,
                error_out = reconstruction_error,
                method = "nnls"
            ))
        }


        ## Multiplicative-update projection

        upper <- max(X)
        lower <- min(X)

        if (!is.finite(upper) || upper <= 0) {
            upper <- 1
        }

        S <- matrix(
            runif(
                ncol(G) * ncol(X),
                min = lower,
                max = upper
            ),
            nrow = ncol(G),
            ncol = ncol(X)
        )

        error_history <- numeric(max_ite)

        for (iteration in seq_len(max_ite)) {
            numerator_S <- crossprod(G, X)

            denominator_S <-
                crossprod(G, G %*% S)

            S <- S * numerator_S / denominator_S
            S[!is.finite(S)] <- 0

            error_history[iteration] <- sum(
                (X - G %*% S)^2
            )

            if (iteration > 1L) {
                previous_error <- error_history[iteration - 1L]
                current_error <- error_history[iteration]

                relative_change <- abs(
                    current_error - previous_error
                ) / max(previous_error, 1e-12)

                if (relative_change <= max_err) {
                    error_history <-
                        error_history[seq_len(iteration)]
                    break
                }
            }
        }

        rownames(S) <- colnames(G)
        colnames(S) <- colnames(X)

        list(
            score = S,
            error_out = error_history,
            method = "multiplicative_update"
        )
    }

    score_list <- vector("list", N)

    ## Projection errors are populated only for datasets that do not
    ## have a singleton second-stage fit.
    projection_error_out <- vector("list", N)
    projection_method <- vector("list", N)

    for (dataset_index in seq_len(N)) {
        score_list[[dataset_index]] <- vector("list", K)

        singleton_steps <- which(vapply(
            group,
            function(grp) {
                length(grp) == 1L &&
                    dataset_index %in% grp
            },
            logical(1)
        ))

        if (length(singleton_steps) > 1L) {
            stop(
                "More than one singleton group was provided for dataset ",
                dataset_index,
                ": group positions ",
                paste(singleton_steps, collapse = ", "),
                "."
            )
        }

        ## Every learned group applicable to this dataset.
        active_groups <- which(vapply(
            group,
            function(grp) {
                dataset_index %in% grp
            },
            logical(1)
        ))

        if (length(active_groups) == 0L) {
            stop(
                "Dataset ",
                dataset_index,
                " does not belong to any group."
            )
        }

        if (length(singleton_steps) == 1L) {
            ## Case 1: scores were learned during singleton fitting

            singleton_step <- singleton_steps[[1]]

            ## At the singleton step, only components already learned
            ## at that point can be present in its S matrix.
            fitted_groups <- active_groups[
                active_groups <= singleton_step
            ]

            S_dataset <- second_stage_S[[singleton_step]]

            expected_rows <- sum(comp_num[fitted_groups])

            if (nrow(S_dataset) != expected_rows) {
                stop(
                    paste0(
                        "Unexpected number of score rows for dataset ",
                        dataset_index,
                        ". Expected ",
                        expected_rows,
                        " but found ",
                        nrow(S_dataset),
                        "."
                    )
                )
            }

            ## With largest-to-smallest ordering, all groups applicable
            ## to the dataset should have been fitted by its singleton.
            missing_groups <- setdiff(active_groups, fitted_groups)

            if (length(missing_groups) > 0L) {
                stop(
                    paste0(
                        "The singleton group for dataset ",
                        dataset_index,
                        " occurs before applicable groups ",
                        paste(missing_groups, collapse = ", "),
                        ". Groups must be ordered from largest to smallest."
                    )
                )
            }

            projection_error_out[[dataset_index]] <- NULL
            projection_method[[dataset_index]] <- "singleton_fit"
        } else {
            ## Case 2: no singleton group
            ##
            ## Project the learned G matrices back onto the original
            ## dataset, holding all G matrices fixed.

            G_dataset <- do.call(
                cbind,
                second_stage_G[active_groups]
            )

            projection_fit <- estimateScoresFixedG(
                X = dataset[[dataset_index]],
                G = G_dataset,
                use_nnls = use_nnls,
                max_ite = max_ite,
                max_err = max_err
            )

            S_dataset <- projection_fit$score

            fitted_groups <- active_groups

            projection_error_out[[dataset_index]] <-
                projection_fit$error_out

            projection_method[[dataset_index]] <-
                projection_fit$method
        }

        ## Split S_dataset into one score matrix per group

        fitted_end_rows <- cumsum(comp_num[fitted_groups])

        fitted_start_rows <- c(
            1L,
            head(fitted_end_rows, -1L) + 1L
        )

        for (group_index in seq_len(K)) {
            fitted_position <- match(
                group_index,
                fitted_groups
            )

            if (!is.na(fitted_position)) {
                rows <- fitted_start_rows[fitted_position]:
                    fitted_end_rows[fitted_position]

                score_list[[dataset_index]][[group_index]] <-
                    S_dataset[rows, , drop = FALSE]
            } else {
                ## Dataset is not part of this group.
                score_list[[dataset_index]][[group_index]] <-
                    matrix(
                        0,
                        nrow = comp_num[group_index],
                        ncol = ncol(dataset[[dataset_index]])
                    )
            }
        }
    }

    score_list <- scoreNameAssign(
        score_list,
        dataset_name,
        group_name
    )

    score_list <- sampleNameAssign(
        score_list,
        sample_name
    )


    ## Give the first-stage output informative names as well

    names(first_stage_G)   <- dataset_name
    names(first_stage_S)   <- dataset_name
    names(first_stage_fit) <- dataset_name

    for (dataset_index in seq_len(N)) {
        rownames(first_stage_G[[dataset_index]]) <- gene_name

        component_names <- paste0(
            dataset_name[dataset_index],
            "_first_stage_component_",
            seq_len(ncol(first_stage_G[[dataset_index]]))
        )

        colnames(first_stage_G[[dataset_index]]) <- component_names
        rownames(first_stage_S[[dataset_index]]) <- component_names
        colnames(first_stage_S[[dataset_index]]) <-
            sample_name[[dataset_index]]
    }


    ## Project an optional additional dataset

    proj_score_list <- NULL
    proj_error_out <- NULL
    proj_method <- NULL

    if (!is.null(proj_dataset)) {
        if (is.null(proj_group)) {
            stop(
                "`proj_group` must be provided when `proj_dataset` is provided."
            )
        }

        if (
            !is.logical(proj_group) ||
            length(proj_group) != length(group) ||
            anyNA(proj_group)
        ) {
            stop(
                "`proj_group` must be a non-missing logical vector of length ",
                length(group),
                "."
            )
        }

        if (!any(proj_group)) {
            stop("At least one element of `proj_group` must be TRUE.")
        }

        proj_dataset <- as.matrix(proj_dataset)

        if (any(!is.finite(proj_dataset))) {
            stop("`proj_dataset` must contain only finite values.")
        }

        if (any(proj_dataset < 0)) {
            stop("`proj_dataset` must be nonnegative.")
        }

        ## Save sample names before preprocessing.
        proj_sample_name <- sampleNameExtractor(proj_dataset)

        ## Work on a separate copy so that feature matching does not alter
        ## the returned linked_component_list.
        projection_component_list <- linked_component_list

        component_gene_names <-
            rownames(projection_component_list[[1]])

        projection_gene_names <- rownames(proj_dataset)

        if (
            !is.null(component_gene_names) &&
            !is.null(projection_gene_names)
        ) {
            ## Preserve the ordering of the learned component matrices.
            common_gene_names <- component_gene_names[
                component_gene_names %in% projection_gene_names
            ]

            if (length(common_gene_names) == 0L) {
                stop(
                    "`proj_dataset` has no named features in common with ",
                    "the learned component matrices."
                )
            }

            proj_rows <- match(
                common_gene_names,
                projection_gene_names
            )

            proj_dataset <- proj_dataset[
                proj_rows,
                ,
                drop = FALSE
            ]

            for (group_index in seq_along(projection_component_list)) {
                component_rows <- match(
                    common_gene_names,
                    rownames(projection_component_list[[group_index]])
                )

                if (anyNA(component_rows)) {
                    stop(
                        "The learned component matrices do not have ",
                        "consistent feature names."
                    )
                }

                projection_component_list[[group_index]] <-
                    projection_component_list[[group_index]][
                        component_rows,
                        ,
                        drop = FALSE
                    ]
            }

            message(
                "Input ",
                length(projection_gene_names),
                " features in `proj_dataset`; found ",
                length(common_gene_names),
                " features in common."
            )
        } else {
            expected_rows <- nrow(projection_component_list[[1]])

            if (nrow(proj_dataset) != expected_rows) {
                stop(
                    "`proj_dataset` has ",
                    nrow(proj_dataset),
                    " rows, but the learned components have ",
                    expected_rows,
                    " rows. Supply matching row names or matrices with ",
                    "identical row order and dimensions."
                )
            }
        }

        ## Apply the same normalization convention used for the input data.
        proj_dataset <- normalizeData(
            list(proj_dataset),
            enable_normalization,
            column_sum_normalization,
            nonnegative_normalization = TRUE
        )[[1]]

        active_proj_groups <- which(proj_group)

        ## Combine only the requested fixed component matrices.
        G_projection <- do.call(
            cbind,
            projection_component_list[active_proj_groups]
        )

        projection_fit <- estimateScoresFixedG(
            X = proj_dataset,
            G = G_projection,
            use_nnls = use_nnls,
            max_ite = max_ite,
            max_err = max_err
        )

        S_projection <- projection_fit$score
        proj_error_out <- projection_fit$error_out
        proj_method <- projection_fit$method

        expected_score_rows <- sum(
            comp_num[active_proj_groups]
        )

        if (nrow(S_projection) != expected_score_rows) {
            stop(
                "Projected score matrix has ",
                nrow(S_projection),
                " rows; expected ",
                expected_score_rows,
                "."
            )
        }


        ## Split the combined score matrix into one matrix per group

        proj_score_list <- vector("list", K)

        active_end_rows <- cumsum(
            comp_num[active_proj_groups]
        )

        active_start_rows <- c(
            1L,
            head(active_end_rows, -1L) + 1L
        )

        for (group_index in seq_along(group)) {
            active_position <- match(
                group_index,
                active_proj_groups
            )

            if (!is.na(active_position)) {
                rows <- active_start_rows[active_position]:
                    active_end_rows[active_position]

                proj_score_list[[group_index]] <-
                    S_projection[rows, , drop = FALSE]
            } else {
                ## Preserve jointNMF's structural-zero convention.
                proj_score_list[[group_index]] <-
                    matrix(
                        0,
                        nrow = comp_num[group_index],
                        ncol = ncol(proj_dataset)
                    )
            }
        }

        ## Use the same naming helpers as jointNMF.
        proj_score_list <- scoreNameAssignProj(
            proj_score_list,
            group_name
        )

        proj_score_list <- sampleNameAssignProj(
            proj_score_list,
            proj_sample_name
        )
    } else if (!is.null(proj_group)) {
        warning(
            "`proj_group` was supplied without `proj_dataset` and will be ignored."
        )
    }

    return(list(
        linked_component_list = linked_component_list,
        score_list = score_list,
        proj_score_list = proj_score_list,
        error_out = second_stage_error,

        first_stage = list(
            component_list = first_stage_G,
            score_list = first_stage_S,
            fit_list = first_stage_fit
        ),

        projection_error_out = projection_error_out,
        projection_method = projection_method,

        proj_error_out = proj_error_out,
        proj_method = proj_method,

        mask_list = second_stage_mask,
        target_info = second_stage_target
    ))
}
