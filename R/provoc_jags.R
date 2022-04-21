#' Convert a list of MCMC objects to a data frame with `chain` and `iter`
#' 
#' @param mcmc.list The list of MCMC objects, as returned by runjags or rjags.
#' @param pivot Use tidyr::pivot_longer to make an easy-to-plot object.
#' 
#' @return A data frame, either with columns iter, chain and all monitored variables if pivot = FALSE or iter, chain, value for the posterior sample and name for the name of the monitored variable if pivot = TRUE.
#' @export
melt_mcmc <- function(mcmc.list, var_names = NULL, pivot = FALSE) {
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
#' @param return_df If TRUE, returns a df with the summary statistics. \code{convergence(res)} will return the convergence information.
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
provoc_jags <- function(
        coco, varmat, 
        adapt = 500, burnin = 1000, thin = 4, 
        sample = 1000, n.chains = 3,
        quiet = TRUE, summarise = FALSE, return_df = TRUE){
    muts <- coco$mutation
    cou2 <- coco$count
    cov2 <- coco$coverage
    vari2 <- varmat

    jags_rho_inits <- function(varmat) {
        method = round(runif(1))
        if(method == 1) {
            init1 <- rho_initializer(varmat)
            init1 <- init1 + rnorm(length(init1), 0, 0.05)
            return(to_feasible(init1))
        } else {
            init2 <- runif(nrow(varmat), 0.2, 0.8)
            return(to_feasible(init2))
        }
    }

    res_temp <- tryCatch(
            runjags::run.jags(
                model = system.file("extdata/provoc.JAGS",
                 package = "provoc"),
                data = list(
                    count = cou2,
                    N = length(cou2),
                    coverage = cov2,
                    P = nrow(vari2),
                    variantmat = vari2,
                    alpha = 2,
                    beta = 8),
                inits = list(
                    list(p1 = jags_rho_inits(varmat), 
                        .RNG.name="base::Super-Duper", 
                        .RNG.seed = 1), 
                    list(p1 = jags_rho_inits(varmat), 
                        .RNG.name="base::Wichmann-Hill", 
                        .RNG.seed = 2), 
                    list(p1 = jags_rho_inits(varmat), 
                        .RNG.name="base::Mersenne-Twister", 
                        .RNG.seed = 3)
                ),
                adapt = adapt,
                burnin = burnin,
                n.chains = n.chains,
                sample = sample,
                thin = thin,
                monitor = c("p"),
                summarise = FALSE,
                silent = quiet,
                #silent.runjags = TRUE,
                method = "parallel"),
            error = function(e) e)
    if("error" %in% class(res_temp)) {
        return(res_df = data.frame(note = res_temp, convergence = FALSE, convergence_note = ""))
    }
    for(i in 1:length(res_temp$mcmc)) {
        colnames(res_temp$mcmc[[i]]) <- rownames(vari2)
    }

    if(!return_df) {
        return(res_temp)
    } else {

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
        convergence_note <- ifelse(convergence_gr,
            yes = g_diag$mpsrf,
            no = paste0("Upper CI larger than 1.15: ", 
                paste(rownames(gr)[gr[, 2] >= 1.15], collapse = ", ")))
        return(list(res_df = res_df, convergence = convergence_gr, convergence_note = convergence_note))
    }
}
