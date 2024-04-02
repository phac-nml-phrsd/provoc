#' Transform a vector to the interior of the feasible region to be used as inits for constrOptim
#' 
#' @param x A vector of numbers
#' 
#' @return A vector with the same length of x such that the sum of the values is larger than 0 but less than 1 and all values are between 0 and 1
to_feasible <- function(x) {
    x <- x - min(x, na.rm = TRUE) + 0.0001
    x <- 0.99 * x / sum(x, na.rm = TRUE)
    x
}


#' Initialize a vector of proportions for the variants of concern, prioritizing current (March 2022) most probable VoCs.
#' 
#' @param varmat The variant matrix to be used in the study; only used for the rownames (which should be VOCs in the expected format)
#' 
#' @return varmat a vector with the same length as the number of rows of varmat, such that the values sum to less than one and each value is between 0 and 1
rho_initializer <- function(varmat) {
    rho_init <- rep(10, length = nrow(varmat))
    # Omicron
    probables1 <- c("BA.1$", "BA.2$", "^B.1.1.529$")
    for (i in seq_along(probables1)) {
        var_in_names <- grepl(probables1[i], rownames(varmat))
        if (any(var_in_names)) {
            rho_init[which(var_in_names)] <- 50
        }
    }

    # Delta
    probables2 <- grepl("AY*", rownames(varmat))
    rho_init[which(probables2)] <- 25

    # Ensure it's not on the boundary
    0.99 * rho_init / sum(rho_init)
}


#' Estimate the proportions of VOCs using Constrained Optimization
#' 
#' If 
#' 
#' @param coco A data frame containing columns labelled count, coverage, and mutation.
#' @param varmat The variant matrix to be used in the study. The rownames must be the VoCs and the colnames must be the mutation names (in the same format as the mutation names in `coco`)
#' @param bootstrap_samples The number of bootstrap samples to use.
#' @param verbose Print messages to the console, default to True.
#' 
#' @return A list containing the results as well as convergence information from \code{constrOptim}.
#' 
#' \describe{
#'      \item{res_df}{The estimated proportions of the variants of concern (\eqn{rho}), including CI if \code{bootstrap_samples > 0}}
#'      \item{convergence}{Logical}
#'      \item{convergence_note}{Convergence code from \code{constrOptim}. Also includes the method used for initializing \eqn{\rho}.}
#' }
#' 
#' @export
#' 
#' @details The estimates are found by minimizing the squared difference between the frequency of each mutation and the prediction of a binomial model where the proportion is equal to the sum of rho times the relevant column of varmat and the size parameter is equal to the coverage. 
#' 
#' The algorithm will first try a prior guess based on the current (March 2022) most common VOCs, then will try a uniform proportion, then (the nuclear option) will try 20 random perturbations until it works. Fails gracefully, with list elements indicating the convergence status and the initialization of rho, and returns the results that had the lowest value of the objective function. 
#' 
#' Bootstrapping is performed parametrically, assuming that coverage is Poisson with a mean of the observed coverage and, based on the sampled value of the coverage, the count is sampled from a binomial distribution with proportion equal to that in the data.
#' 
#' @seealso \code{\link[stats]{constrOptim}}
#' 
#' @examples
#' varmat <- simulate_varmat() # default values (Omicron)
#' varmat <- varmat[row.names(varmat) %in% c("B.1.1.529", "BA.1", "BA.2")]
#' coco <- simulate_coco(varmat, rel_counts = c(100, 200, 300)) # expect 1/6, 2/6, and 3/6
#' res <- copt_binom(coco, varmat)
#' res$res_df
provoc_optim <- function(coco, varmat, bootstrap_samples = 0,
    verbose = TRUE, rho_init = NULL) {
    muts <- coco$mutation[coco$coverage > 0]
    cou2 <- coco$count[coco$coverage > 0]
    cov2 <- coco$coverage[coco$coverage > 0]
    vari2 <- varmat[, coco$coverage > 0]
    if (is.null(rho_init)) {
        rho_init <- rho_initializer(vari2)
    } else {
        rho_init <- to_feasible(rho_init)
    }

    objective <- function(rho, count, varmat, coverage) {
        -sum(stats::dbinom(x = as.numeric(count),
                size = as.numeric(coverage),
                prob = 0.999 * rho %*% varmat + 0.0001,
                log = TRUE))
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    ui <- rbind(u_sum1, u_p0)
    ci <- c(c_sum1, c_p0)


    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        count = cou2, coverage = cov2, varmat = vari2,
        control = list(maxit = 10000))
    res$init_method <- "Prior_Assumption"

    if (res$convergence) { # if not converged,
        # Try a different initialization
        rho_init <- 0.99 * rep(1 / nrow(vari2), nrow(vari2))
        res <- stats::constrOptim(rho_init,
            f = objective, grad = NULL,
            ui = ui, ci = ci,
            count = cou2, coverage = cov2, varmat = vari2,
            control = list(maxit = 10000))
        res$init_method <- "Uniform"
    }

    bestres <- res
    if (res$convergence) {
        if (verbose) {
            print("Trying the nuclear option for constrOptim.")
            print("This has never actually worked before.")
            print("Godspeed.")
        }
        # Uniform inititialization
        i <- 0
        converged <- FALSE
        while (i < 20 && !converged) {
            i <- i + 1
            if (!i %% 10) print(paste0("Attempt ", i, " of 20."))

            # Add noise to previous iteration
            rho_init <- res$par +
                stats::rnorm(length(rho_init), 0, 0.05)
            # Constrain to interior of feasible region
            rho_init <- to_feasible(rho_init)

            res <- stats::constrOptim(rho_init,
                f = objective, grad = NULL,
                ui = ui, ci = ci,
                count = cou2, coverage = cov2, varmat = vari2,
                control = list(maxit = 1000))
        }
        res$init_method <- "Nuclear"

        # Only take the best results out of all nuclear tries
        if (res$value < bestres$value) {
            bestres <- res
        }

        converged <- !res$convergence
        if (i == 20 && verbose) {
            print("Nuclear Option failed; going with best results.")
        }
    }

    res_df <- data.frame(rho = bestres$par,
        ci_low = NA,
        ci_high = NA,
        variant = rownames(varmat))
    convergence <- ifelse(bestres$convergence == 0, TRUE, bestres$convergence)
    convergence_note <- paste("Optim results: ",
        bestres$convergence,
        "; Initialization: ", bestres$init_method, sep = "")
    boots <- NULL
    if (bootstrap_samples > 0) {
        resampled_coverage <- rmultinom(bootstrap_samples,
            size = sum(cov2),
            prob = cov2 / sum(cov2)) |>
            as.integer()
        resampled_counts <- rbinom(length(resampled_coverage),
            size = resampled_coverage,
            prob = rep(cou2 / cov2, bootstrap_samples))
        resamples <- data.frame(
            count = resampled_counts,
            coverage = resampled_coverage,
            mut = rep(muts, bootstrap_samples),
            iteration = rep(1:bootstrap_samples,
                each = length(muts))
        )

        boots <- sapply(split(resamples, resamples[, "iteration"]),
            function(x) {
                provoc_optim(x, vari2,
                    rho_init = res_df$rho)$res_df[, "rho"]
            }
        ) |> t()
        colnames(boots) <- res_df$variant

        ci <- apply(boots, 2, quantile, prob = c(0.025, 0.975))
        res_df$ci_low <- ci[1, ]
        res_df$ci_high <- ci[2, ]
    }

    return(list(res_df = res_df,
            convergence = convergence,
            convergence_note = convergence_note,
            bootstrap_samples = boots))
}
