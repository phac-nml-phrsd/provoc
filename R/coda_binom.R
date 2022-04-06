#' Convert a list of MCMC objects to a data frame with `chain` and `iter`
#' 
#' @param mcmc.list The list of MCMC objects, as returned by runjags or rjags
#' @param pivot Use tidyr::pivot_longer to make an easy-to-plot object
#' 
#' @return A data frame, either with columns iter, chain and all monitored variables if pivot = FALSE or iter, chain, value for the posterior sample and name for the name of the monitored variable if pivot = TRUE.
#' @export
melt_mcmc <- function(mcmc.list, pivot = TRUE) {
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

#' Estimate a Bayesian model for rho.
#' 
#' @param coco A data frame containing counts, coverage, and mutation names
#' @param varmat The variant matrix to be used in the study. The rownames must be the VoCs and the colnames must be the mutation names (in the same format as the mutation names in `coco`)
#' @param prm The parameters for the mcmc fit, generally from a prm.json file.
#' 
#' @return an mcmc.list object with each column representing the proportion of the variant of concern, in the order of the rownames of variantmat. It is suggested to use melt_mcmc to get output that plays better with ggplot2 and dplyr.
#' @export
coda_binom <- function(
        coco, 
        varmat, 
        prm = list(adapt = 500, burnin = 1000, 
            sample = 1000, thin = 4, quiet = 0)){
    if(requireNamespace("runjags", quietly = TRUE)) {
        res <- tryCatch(
                runjags::run.jags(
                    model = system.file("inst/extdata/provoc.JAGS", package = "provoc"),
                    data = list(
                        count = coco$count,
                        N = nrow(coco),
                        coverage = coco$coverage + 1,
                        P = nrow(varmat),
                        variantmat = varmat,
                        alpha = 2,
                        beta = 8),
                    #inits = list(p1 = rep(1/nrow(variantmat), nrow(variantmat))),
                    adapt = prm$adapt,
                    burnin = prm$burnin,
                    n.chains = 3,
                    sample = prm$sample,
                    thin = prm$thin,
                    monitor = c("p"),
                    summarise = FALSE,
                    silent = prm$quiet,
                    #silent.runjags = TRUE,
                    method = "parallel"),
                error = function(e) e)
        return(res)
    } else {
        "JAGS and runjags are required to use this function."
    }
}
