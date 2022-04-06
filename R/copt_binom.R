#' Transform a vector to the interior of the feasible region to be used as inits for constrOptim
#' 
#' @param x A vector of numbers
#' 
#' @return A vector with the same length of x such that the sum of the values is larger than 0 but less than 1 and all values are between 0 and 1
to_feasible <- function(x){
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
    # B or Omicron
    probables1 <- c("BA.1", "BA.2", "B.1.1.529", "B")
    for (i in seq_along(probables1)) {
        if(probables1[i] %in% rownames(varmat)) {
            row_index <- which(
                rownames(varmat) == probables1[i]
            )
            rho_init[row_index] <- 50
        }
    }

    # Delta
    probables2 <- grepl("^AY", rownames(varmat))
    rho_init[probables2] <- 25

    # Ensure it's not on the boundary
    0.99 * rho_init / sum(rho_init)
}


#' Estimate the proportions of VOCs using Constrained OPTimization
#' 
#' @param coco A data frame containing counts, coverage, and mutation names
#' @param varmat The variant matrix to be used in the study. The rownames must be the VoCs and the colnames must be the mutation names (in the same format as the mutation names in `coco`)
#' 
#' @return A list containing the vector of proportions for the variants of concern, in the same order as the rows of varmat. This list also contains valuable convergence information.
#' @export
#' 
#' @details The estimates are found by minimizing the squared difference between the frequency of each mutation and the prediction of a binomial model where the proportion is equal to the sum of rho times the relevant column of varmat and the size parameter is equal to the coverage. 
#' 
#' The algorithm will first try a prior guess based on the current (March 2022) most common VOCs, then will try a uniform proportion, then (the nuclear option) will try 100 random perturbations until it works. Fails gracefully, with list elements indicating the convergence status and the initialization of rho. 
#' 
#' This function currently does not return any estimate of variance for the proportions and should not be trusted beyond a quick check.
copt_binom <- function(coco, varmat) {
    bad_freq <- which(is.na(coco$frequency))
    muts <- coco$mut[-bad_freq]
    freq2 <- coco$frequency[-bad_freq]
    cov2 <- coco$coverage[-bad_freq]
    vari2 <- varmat[, muts]

    objective <- function(rho, frequency, varmat, coverage) {
        count <- round(frequency * coverage)
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        -sum(stats::dbinom(x = count, size = coverage, prob = prob,
            log = TRUE))
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # sum(p) > 0
    u_sum0 <- rep(1, length(rho_init))
    c_sum0 <- 0

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    # p_i < 1 => -p_i > -1 => -p_i + 1 > 0
    u_p1 <- diag(length(rho_init))
    c_p1 <- rep(-1, length(rho_init))

    ui <- rbind(u_sum1, u_sum0, u_p0, u_p1)
    ci <- c(c_sum1, c_sum0, c_p0, c_p1)


    rho_init <- rho_initializer(vari2)
    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        frequency = freq2, coverage = cov2, varmat = vari2,
        control = list(maxit = 100000))
    res$init_method <- "Prior_Assumption"

    if(res$convergence) { # if not converged,
        # Try a different initialization
        rho_init <- 0.99*rep(1/nrow(vari2), nrow(vari2))
        res <- stats::constrOptim(rho_init,
            f = objective, grad = NULL,
            ui = ui, ci = ci,
            frequency = freq2, coverage = cov2, varmat = vari2,
            control = list(maxit = 100000))
        res$init_method <- "Uniform"
    }

    if(res$convergence) { 
        print("Trying the nuclear option for constrOptim.")
        print("Godspeed.")
        # Uniform inititialization
        i <- 0
        converged <- FALSE
        while(i < 100 & !converged) {
            i <- i + 1
            print(paste0("Attempt ", i, " of 100."))
            # Add noise to previous iteration
            rho_init <- res$par +
                stats::rnorm(length(rho_init), 0, 0.1)
            # Constrain to interior of feasible region
            rho_init <- to_feasible(rho_init)

            res <- stats::constrOptim(rho_init,
                f = objective, grad = NULL,
                ui = ui, ci = ci,
                frequency = freq2, coverage = cov2, varmat = vari2,
                control = list(maxit = 10000))
            }
            res$init_method <- "Nuclear"

            converged <- !res$convergence
    }

    res
}
