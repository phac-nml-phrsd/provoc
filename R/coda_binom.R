#' Convert a list of MCMC objects to a data frame with `chain` and `iter`
#' 
#' @param mcmc.list The list of MCMC objects, as returned by runjags or rjags.
#' @param pivot Use tidyr::pivot_longer to make an easy-to-plot object.
#' 
#' @return A data frame, either with columns iter, chain and all monitored variables if pivot = FALSE or iter, chain, value for the posterior sample and name for the name of the monitored variable if pivot = TRUE.
#' @export
melt_mcmc <- function(mcmc.list, pivot = FALSE) {
    mcmc <- dplyr::bind_rows(lapply(1:length(mcmc.list$mcmc),
            function(i) {
                x <- as.data.frame(mcmc.list$mcmc[[i]])
                x$chain <- i
                x$iter <- 1:nrow(x)
                x
            }
        ))
    if(pivot) {
        return(tidyr::pivot_longer(mcmc, -c("chain", "iter")))
    } else {
        return(mcmc)
    }
}

#' Estimate a Bayesian model for \eqn{\rho}.
#' 
#' Fits \eqn{counts_i \sim Binom(p_i, coverage_i)}, where \eqn{p_i} is \eqn{\rho} (the vector of proportions for each variant) times the column of the variant matrix corresponding to the \eqn{i}th mutation. In other words, \eqn{p_i} is the sum of the proportions of the variants that contain mutation i.
#' 
#' @param coco A data frame containing counts, coverage, and mutation names
#' @param varmat The variant matrix to be used in the study. The rownames must be the VoCs and the colnames must be the mutation names (in the same format as the mutation names in `coco`)
#' @param adapt,burnin,sample,thin Parameters passed to runjags. Note that \code{sample} is the final number of samples that you will receive from each chain (i.e. if \code{thin = 5} and \code{sample = 1000}, you will end up drawing 5000 samples so that 1000 are returned after thinning).
#' @param quiet If TRUE, much (but not all) of the runjags output will be suppressed.
#' 
#' @return an mcmc.list object with each column representing the proportion of the variant of concern, in the order of the rownames of variantmat. It is suggested to use melt_mcmc to get output that plays better with ggplot2 and dplyr.
#' @export
#' 
#' @examples
#' varmat <= simulate_varmat()
#' coco <- simulate_coco(varmat, rel_counts = c(100, 200, 300))
#' 
#' coda <- coda_binom(coco, varmat)
#' codf <- melt_mcmc(coda, pivot = FALSE)
#' 
#' plot(x = codf$iter, y = codf[, 1], col = codf$chain)
#' 
coda_binom <- function(
        coco, varmat, 
        adapt = 500, burnin = 1000, sample = 1000, thin = 4,
        quiet = TRUE){
    if(requireNamespace("runjags", quietly = TRUE)) {
        bad_freq <- which(is.na(coco$coverage))
        if(length(bad_freq) > 0) {
            muts <- coco$mutation[-bad_freq]
            cou2 <- coco$count[-bad_freq]
            cov2 <- coco$coverage[-bad_freq]
            vari2 <- varmat[, muts]
        } else {
            muts <- coco$mutation
            cou2 <- coco$count
            cov2 <- coco$coverage
            vari2 <- varmat
        }
        res <- tryCatch(
                runjags::run.jags(
                    model = system.file("extdata/provoc.JAGS", package = "provoc"),
                    data = list(
                        count = cou2,
                        N = length(cou2),
                        coverage = cov2 + 1,
                        P = nrow(vari2),
                        variantmat = vari2,
                        alpha = 2,
                        beta = 8),
                    #inits = list(p1 = rep(1/nrow(variantmat), nrow(variantmat))),
                    adapt = adapt,
                    burnin = burnin,
                    n.chains = 3,
                    sample = sample,
                    thin = thin,
                    monitor = c("p"),
                    summarise = FALSE,
                    silent = quiet,
                    #silent.runjags = TRUE,
                    method = "parallel"),
                error = function(e) e)
        return(res)
    } else {
        "JAGS and runjags are required to use this function."
    }
}
