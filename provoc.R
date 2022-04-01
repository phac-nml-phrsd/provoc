# Functions for ProVoC

p_model <- {"
model{
    # Likelihood
    for(i in 1:N) {
        # Ensure that the proportion is less than 1
        #p1[i] <- min(0.9999, pvec[i])
        #sump[i] <- sum(pvec[i])
        count[i] ~ dbinom((pvec[i] + 0.001)*0.999, coverage[i])
    }

    # Prior on proportions for each variant
    for(j in 1:P) {
        p1[j] ~ dbeta(alpha, beta)
        p[j] <- ifelse(sum(p1) >= 1, p1[j]/sum(p1), p1[j])
    }

    # Sum of p times an indicator function
    pvec <- p %*% variantmat
}
"}



melt_mcmc <- function(mcmc.list, pivot = TRUE) {
    mcmc <- bind_rows(lapply(1:length(mcmc.list$mcmc),
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



to_feasible <- function(x){
    x <- x - min(x, na.rm = TRUE) + 0.0001
    x <- 0.99 * x / sum(x, na.rm = TRUE)
    x
}



coda_binom <- function(cocovoc1, vari2, prm){
    tryCatch(
            run.jags(
                model = p_model,
                data = list(
                    count = cocovoc1$count,
                    N = nrow(cocovoc1),
                    coverage = cocovoc1$coverage + 1,
                    P = nrow(vari2),
                    variantmat = vari2,
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
}



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



copt_binom <- function(coco, varmat) {
    bad_freq <- which(is.na(coco$frequency))
    muts <- coco$mut[-bad_freq]
    freq2 <- coco$frequency[-bad_freq]
    cov2 <- coco$coverage[-bad_freq]
    vari2 <- varmat[, muts]
    rho_init <- 0.99*rep(1/nrow(vari2), nrow(vari2))

    objective <- function(rho, frequency, varmat, coverage) {
        count <- round(frequency * coverage)
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        -sum(dbinom(x = count, size = coverage, prob = prob,
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
    res <- constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        frequency = freq2, coverage = cov2, varmat = vari2,
        control = list(maxit = 100000))
    res$init_method <- "Prior_Assumption"

    if(res$convergence) { # if not converged,
        # Try a different initialization
        res <- constrOptim(rho_init,
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
                rnorm(length(rho_init), 0, 0.1)
            # Constrain to interior of feasible region
            rho_init <- to_feasible(rho_init)

            res <- constrOptim(rho_init,
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


