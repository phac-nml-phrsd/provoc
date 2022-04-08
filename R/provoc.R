#' Proportions of Variants of Concern
#' 
#' Un-fuses coco and varmat and applies the appropriate estimation technique. Applies to each unique value of 
#' 
#' @param fused The fused data frame of coco and varmat
#' @param method Optimization or Bayesian?
#' @param ... Arguments to be passed to provoc_jags (ignored for method = "optim"). See \code{?provoc_binom}.
#' 
#' @return Results of the estimation. Currently has different output depending on method, but will hopefully be unified in a later version.
#' @export
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
    names(res_list) <- samples
    for(i in seq_along(res_list)) {
        # TODO: Allow parameters to be passed to the respective functions.
        cat("\n")
        message(paste0("Fitting sample ", samples[i], ", ", which(samples == samples[i]), "of ", length(samples)))
        fusi <- fused[fused$sample == samples[i], ]
        variants <- startsWith(names(fusi), "var_")
        varnames <- names(fused)[variants]
        coco <- fusi[, !variants]

        vardf <- fusi[, variants]
        varmat <- t(as.matrix(vardf))
        varmat <- matrix(as.numeric(varmat), ncol = ncol(varmat))
        rownames(varmat) <- gsub("var_", "", varnames)
        colnames(varmat) <- coco$mutation

        if(method[1] == "optim") {
            res_temp <- provoc_optim(coco, varmat)
            res <- list(
                point_est = data.frame(rho = res_temp$par, 
                    variant = varnames,
                    ci_low = NA,
                    ci_high = NA),
                convergence = res_temp$convergence,
                convergence_note = res_temp$init_method,
                note = "CI is NA; bootstrapping not yet implemented.",
                logLik = -res_temp$value)
        } else {
            if(requireNamespace("runjags", quietly = TRUE)) {
                res_temp <- provoc_jags(coco, varmat, ...)

                if(requireNamespace("coda", quietly = TRUE)) {
                    gr <- coda::gelman.diag(res_temp)$psrf
                    convergence_gr <- coda::gelman.diag(res_temp)$mpsrf > 1.15
                } else {
                    gr <- "coda not install, cannot calculate GR statistics"
                    convergence_gr <- NA
                }

                point_est <- melt_mcmc(res_temp, varmat, pivot = FALSE)
                point_est <- point_est[, -which(names(point_est) %in% c("chain", "iter"))]
                point_est <- t(apply(point_est, 2, quantile, probs = c(0.025, 0.5, 0.975)))
                point_est <- as.data.frame(point_est)
                names(point_est) <- c("ci_low", "rho", "ci_high")
                point_est$variant <- rownames(varmat)
                point_est <- point_est[, c("rho", "ci_low", "ci_high", "variant")]
                rownames(point_est) <- NULL

                logLik <- sum(stats::dbinom(
                    x = fused$count, 
                    size = fused$coverage, 
                    prob = as.numeric(point_est$rho %*% varmat), 
                    log = TRUE))

                res <- list(
                    point_est = point_est,
                    convergence = convergence_gr,
                    convergence_note = gr,
                    note = "To return entire posteriors, call provoc_binom() directly. Loglik is calculated for medians only. Convergence is based on multivariate psrf.",
                    logLik = logLik)
                    
            } else {
                res <- "JAGS and runjags are required to use this method."
            }
        }
        # TODO: (Long Term) make this it's own class with nice printing defaults
        res_list[[i]] <- res

        message("Done.")
    }

    res_list
}



