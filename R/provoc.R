#' Proportions of Variants of Concern
#' 
#' Un-fuses coco and varmat and applies the appropriate estimation technique. If a column labelled "sample" is present, applies the analysis to each sample separately.
#' 
#' @param fused The fused data frame of coco and varmat. The fusion ensures that the mutations are properly joined and varmat contains only the most pertinent mutations.
#' @param method Optimization or Bayesian?
#' @param ... Arguments to be passed to \code{provoc_jags()} (ignored for method = "optim"). See \code{?provoc_jags}.
#' 
#' @return Results of the estimation. Regardless of the technique used, the results are as follows.
#' 
#' \describe{
#'      \item{point_est}{A data frame with columns labelled rho, ci_low, ci_high, and the name of the variant. For \code{optim}, these are the estimates (bootstrap CI not yet implemented). For \code{runjags}, these are the medians and quantiles.}
#'      \item{convergence}{True if the algorithm converged.}
#'      \item{convergence_notes}{For \code{optim}, the message returned from \code{optim} as well as the initialization method for rho. For \code{runjags}, the gelman-rubin statistics for each parameter.}
#'      \item{note}{A brief note about the method.}
#'      \item{loglik}{The log likelihood value. For \code{runjags}, this is calculated for the medians of all parameters. This may be changed in future updates.}
#'      \item{sample_info}{The function checks for any columns that have exactly one unique value within a given sample. The unique value of each is returned in a data frame. Very useful if samples have columns such as "date" or "location".}
#' }
#' @export
#' 
#' @examples
#' varmat <- simulate_varmat()
#' coco <- simulate_coco(varmat)
#' fused <- fuse(coco, varmat)
#' res <- provoc(fused)
#' res$point_est
provoc <- function (fused, method = c("optim", "runjags"), ...) {
    if("sample" %in% colnames(fused)) {
        samples <- unique(fused$sample)
        sample_table <- table(fused$sample)
        if(any(sample_table < 5)) {
            if(mean(sample_table < 5) < 0.5) {
                print(sample_table)
                stop("Too few samples")
            } else {
                samples <- names(sample_table[sample_table > 5])
                warning("Some samples have fewer than 5 observations and have been removed from analysis.")
                print(table(samples))
            }
        }
        if(any(sample_table < 10)) {
            warning("At least one of the samples has fewer than 10 observations.")
            print(sample_table)
        }
    } else {
        fused$sample <- 1
        samples <- 1
    }

    res_list <- vector(mode = "list", length = length(samples))
    convergence_note <- character(length(samples))
    convergence <- logical(length(samples))
    names(res_list) <- samples
    for(i in seq_along(res_list)) { # TODO: Parallelize (for optim)
        cat("\n")
        message(paste0("Fitting sample ", samples[i], ", ", 
            which(samples == samples[i]), " of ", length(samples)))
        fusi <- fused[fused$sample == samples[i], ]
        variants <- startsWith(names(fusi), "var_")
        varnames <- names(fused)[variants]
        coco <- fusi[, !variants]

        # Guaranteed to include a column called "sample"
        sample_info <- apply(coco, 2, function(x) length(unique(x)) == 1)
        sample_info <- coco[1, sample_info, drop = FALSE]

        vardf <- fusi[, variants]
        varmat <- t(as.matrix(vardf))
        varmat <- matrix(as.numeric(varmat), ncol = ncol(varmat))
        rownames(varmat) <- gsub("var_", "", varnames)
        colnames(varmat) <- coco$mutation

        if(method[1] == "optim") {
            res_temp <- provoc_optim(coco, varmat)

            res_df <- data.frame(rho = res_temp$par,
                    ci_low = NA,
                    ci_high = NA, 
                    variant = rownames(varmat))
            for (ii in seq_len(ncol(sample_info))) {
                res_df[, names(sample_info)[ii]] <- sample_info[1, ii]
            }
            convergence[i] <- res_temp$convergence == 0
            convergence_note[i] <- paste("Optim results: ", 
                res_temp$convergence, 
                "; Initialization: ", res_temp$init_method, sep = '')
            res_list[[i]] <- res_df
        } else {
            if(requireNamespace("runjags", quietly = TRUE)) {
                res_temp <- provoc_jags(coco, varmat, ...)

                # TODO: all of this should be in the provoc_jags() function, not cluttering up this file.

                if(requireNamespace("coda", quietly = TRUE)) {
                    g_diag <- coda::gelman.diag(res_temp)
                    gr <- g_diag$psrf
                    convergence_gr <- g_diag$mpsrf <= 1.15
                } else {
                    gr <- "coda library not installed, cannot calculate GR statistics"
                    convergence_gr <- NA
                }

                point_est <- melt_mcmc(res_temp, varmat,
                    pivot = FALSE)
                point_est <- point_est[, 
                    -which(names(point_est) %in% c("chain", "iter"))]
                point_est <- t(apply(point_est, 2, quantile, 
                    probs = c(0.025, 0.5, 0.975)))
                point_est <- as.data.frame(point_est)
                names(point_est) <- c("ci_low", "rho", 
                    "ci_high")
                point_est$variant <- rownames(varmat)
                point_est <- point_est[, c("rho", "ci_low", 
                    "ci_high", "variant")]
                rownames(point_est) <- NULL

                res_df <- point_est
                for (ii in seq_len(ncol(sample_info))) {
                    res_df[, names(sample_info)[ii]] <- sample_info[1, ii]
                }
                convergence[i] <- convergence_gr
                convergence_note[i] <- ifelse(convergence_gr,
                    yes = g_diag$mpsrf,
                    no = paste0("Upper CI larger than 1.15: ", 
                        paste(rownames(gr)[gr[, 2] >= 1.15], collapse = ", ")))
                res_list[[i]] <- res_df
                    
            } else {
                res_df <- "JAGS and runjags are required to use this method."
                res_list[[i]] <- res_df
            }
        }
        # TODO: Add methods (for both single results and lists): summary, plot
            # Only print top variants; include columns that are constant within samples; convergence status; log-Likelihood
            # Autoplot (gg) for res and res_list objects?
                # Bonus: colour schemes that respect variant names? Could be another repo that other labs might enjoy.
        # TODO: Replace processing in provoc() with separate processing functions for optim and jags; user can choose how to process.
        # TODO: Include median read depth, quality measures

        message("Done.")
    }

    res <- dplyr::bind_rows(res_list)
    attr(res, "convergence") <- data.frame(sample = samples, 
        convergence = convergence, 
        convergence_note = convergence_note)
    res
}

#' Check if provoc converged
#' 
#' If converged, returns True and prints a message. Otherwise, prints the samples and the note giving hints as to why it didn't converge.
#' 
#' @param res The result of \code{provoc()}
#' @param verbose Print a message to the screen? 
#' 
#' @return Invisbly returns TRUE if all samples converged, false otherwise. 
#' @export
convergence <- function(res, verbose = TRUE) {
    if(!"convergence" %in% attributes(attributes(res))$names) {
        stop("Not a result of provoc - does not have correct attributes")
    }
    conv <- attr(res, "convergence")
    if(any(!as.logical(conv$convergence))) {
        if(verbose) print(conv[which(!as.logical(conv$convergence)), -2])
        return(invisible(FALSE))
    } else {
        if(verbose) cat("All samples converged\n")
        return(invisible(TRUE))
    }
}
