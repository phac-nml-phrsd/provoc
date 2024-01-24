#' Proportions of Variants of Concern
#' 
#' Un-fuses coco and varmat and applies the appropriate estimation technique. If a column labelled "sample" is present, applies the analysis to each sample separately.
#' 
#' @param coco The counts and coverage data frame. Ignored if \code{fused} is supplied.
#' @param varmat The matrix of variant mutations. Ignored if \code{fused} is supplied.
#' @param fused The fused data frame of coco and varmat. The fusion ensures that the mutations are properly joined and varmat contains only the most pertinent mutations.
#' @param ncores For optim, the number of cores to be used. Requires the parallel package for ncores > 1
#' @param update_interval Print message after \code{update_interval} samples. Set to 0 to suppress output.
#' 
#' @return Results of the estimation. Regardless of the technique used, the results are as follows.
#' 
#' \describe{
#'      \item{point_est}{A data frame with columns labelled rho, ci_low, ci_high, and the name of the variant. These are the estimates (bootstrap CI not yet implemented).}
#'      \item{convergence}{True if the algorithm converged.}
#'      \item{convergence_notes}{the message returned from \code{optim} as well as the initialization method for rho.}
#'      \item{note}{A brief note about the method.}
#'      \item{loglik}{The log likelihood value. This may be changed in future updates.}
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
provoc <- function (coco, varmat, fused = NULL, ncores = 1, bootstrap_samples = 0, update_interval = 20, verbose = TRUE) {
    if(is.null(fused)) fused <- fuse(coco, varmat)
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
    t0s <- double(length(samples))
    names(res_list) <- samples
    if (update_interval) message(paste0("Fitting ", length(samples), " samples."))
    for (i in seq_along(res_list)) { # TODO: Parallelize (for optim)
        t0 <- Sys.time()
        if (update_interval){
            if (!i %% update_interval) {
                cat("\n")
                message(paste0("Fitting sample ", samples[i], ", ",
                    which(samples == samples[i]), " of ",
                    length(samples)))
            }
        }
        fissed <- fission(fused, sample = samples[i])
        coco <- fissed$coco
        varmat <- fissed$varmat

        # Guaranteed to include a column called "sample"
        sample_info <- apply(coco, 2, function(x) length(unique(x)) == 1)
        sample_info <- coco[1, sample_info, drop = FALSE]
        
        res_temp <- provoc_optim(coco, varmat, bootstrap_samples = bootstrap_samples, verbose = verbose)
        res_df <- res_temp$res_df
        for (ii in seq_len(ncol(sample_info))) {
            res_df[, names(sample_info)[ii]] <- sample_info[1, ii]
        }
        convergence[i] <- res_temp$convergence
        convergence_note[i] <- res_temp$convergence_note
        res_list[[i]] <- res_df
        # TODO: Add methods (for both single results and lists): summary, plot
            # Only print top variants; include columns that are constant within samples; convergence status; log-Likelihood
            # Autoplot (gg) for res and res_list objects?
                # Bonus: colour schemes that respect variant names? Could be another repo that other labs might enjoy.
        # TODO: Replace processing in provoc() with separate processing functions for optim and jags; user can choose how to process.
        # TODO: Include median read depth, quality measures

        t0s[i] <- difftime(Sys.time(), t0, units = "mins")
    }

    res <- dplyr::bind_rows(res_list)
    attr(res, "convergence") <- data.frame(sample = samples, 
        convergence = convergence, 
        convergence_note = convergence_note,
        time = t0s)
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


#' Predict using Proportions of Variants of Concern
#'
#' Takes a named list with an estimate of the proportions and the
#' associated variant matrix, performs matrix multiplication to
#' predict outcomes, and returns results in the same order as the original data.
#'
#' @param provoc_obj Named list with `proportions` and `variant_matrix`.
#' @param newdata Not yet implemented.
#' @param type Not yet implemented.
#' @param se.fit Not yet implemented.
#' @param dispersion Not yet implemented.
#' @param terms Not yet implemented.
#' @param na.action Not yet implemented.
#' @return Predicted values in the same order as the input data.
#' @export
#' @examples
#' predicted_results <- predict.provoc(provoc_obj)
predict.provoc <- function(provoc_obj, newdata = NULL, type = NULL, se.fit = NULL,
                    dispersion = NULL, terms = NULL, na.action = NULL) {
                        
    # Ensure that provoc_obj is a list and not NULL
    if (is.null(provoc_obj) || !is.list(provoc_obj)) {
        stop("Input 'provoc_obj' must be a non-null list.")
    }
    
    # Check if the required elements are present in the input
    if (!all(c("proportions", "variant_matrix") %in% names(provoc_obj))) {
        stop("Input 'provoc_obj' must contain 'proportions' and 'variant_matrix'.")
    }

    proportions <- provoc_obj$proportions
    variant_matrix <- provoc_obj$variant_matrix

    # Check for numeric types and non-empty inputs
    if (!is.numeric(proportions) || length(proportions) == 0) {
        stop("'proportions' must be a non-empty numeric vector.")
    }
    if (!is.matrix(variant_matrix) || any(dim(variant_matrix) == 0)) {
        stop("'variant_matrix' must be a non-empty numeric matrix.")
    }

    # Error checking for dimensions
    if (ncol(variant_matrix) != length(proportions)) {
        stop("Dimensions of 'variant_matrix' and 'proportions' do not match.")
    }
    
    # Perform matrix multiplication for predictions
    results <- variant_matrix %*% proportions
    
    return(results)  # Return results in the same order
}
