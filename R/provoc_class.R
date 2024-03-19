#' Print the results of lineage abundance estimation
#' @export
print.provoc <- function(provoc_obj, n = 6) {
    cat("Call: ", as.character(attributes(provoc_obj)$formula))
    cat("\n\n")
    cat("Convergence:\n")
    print(unlist(attributes(provoc_obj)$convergence))
    cat("\n")
    n <- min(n, nrow(provoc_obj))
    provoc_df <- provoc_obj[order(-provoc_obj$rho), ]
    provoc_df$rho <- ifelse(provoc_df$rho < 0.001,
        "<0.001", round(provoc_df$rho, 3))
    cat("Top", n, "variants:\n")
    print.data.frame(provoc_df[1:n, ], digits = 3)
}

#' Summarise results of model fitting
#' 
#' Prints the most useful diagnostics to the screen, invisibly returning them as a list.
#' 
#' @export
summary.provoc <- function(provoc_obj) {
    cat("Call: ")
    print(attributes(provoc_obj)$formula)
    cat("\n")
    fitted <- predict(provoc_obj)

    # TODO: Data summary
    # Number of mutations used in the fitting,
    # Similarity of variants (Jaccard, based on data used in model)
    # Entropy of frequencies, or deviation from 0.5
    # Five-number summary for coverage
    # ...

    # TODO: Per-variant goodness of fit diagnostics
    # | variant | total deviance | n_obs | ...
    # |---------|----------------|-------| ...
    # | Overall | ...
    # | B.1.1.7 | ...
    # | B.1.427 | ...

    # TODO: Overall goodness of fit diagnostics
    # R^2 analogue, Chi-squaretest for deviance, number of parameters
    # convergence diagnostics, ...

    # TODO: return a summary.provoc object with associated print method.
    # This is how summary.lm prints to the screen if not assigned to
    # an object, and does not print if it is assigned to an object.
    # This f'n to be split into one function to create the object
    # and a print method.
    return(invisible(list(formula = formula, resids = fitted)))
}

#' Print the summary of a provoc object
print.summary.provoc <- function(summary.provoc) {
    # TODO: All of it.
}

#' Plot the results of model fitting
#' @export
plot.provoc <- function(provoc_obj, plot_type = c("barplot")) {
    #plot_types <- c("b", "r")
    #if(!all(plot_type %in% plot_types)) stop("Invalid plot choice.")

    mfrow <- switch(as.character(length(plot_type)),
        "1" = c(1,1),
        "2" = c(1,2),
        c(2,2))
    par(mfrow = mfrow)

    # Barplot
    if(1 %in% plot_type || any(startsWith(plot_type, "b"))) {
        barplot(
            height = matrix(provoc_obj$rho, 
                ncol = ifelse("group" %in% names(provoc_obj), 
                    length(unique(provoc_obj$group)), 
                    1)), 
            col = 1:length(unique(provoc_obj$variant)),
            horiz = TRUE, 
            xlim = c(0, 1),
            xlab = "Proportion", ylab = NULL)
        legend("topright", 
            legend = unique(provoc_obj$variant), 
            col = 1:length(unique(provoc_obj$variant)), 
            pch = 15)
    }
}

#' Plot the results using ggplot2
#' 
#' Only plots the results, does not plot residuals
#' 
#' @export
autoplot.provoc <- function(provoc_obj, date_col = NULL) {
    if(!"ggplot2" %in% .packages()) {
        stop("Please load ggplot2 before using this function.")
    }

    gg <- ggplot(provoc_obj) + 
            geom_bar(stat = "identity", position = "stack") +
            coord_flip() +
            lims(y = c(0, 1))
    if (!is.null(date_col)) {
        if (!is.Date(provoc_obj[, date_col])) 
            stop("Supplied date column does not include Date values. \nTry lubridate::ymd().")

        gg <- gg  +
            aes(x = 1, y = rho, fill = variant) +
            labs(x = "Proportion", y = NULL, fill = "Lineage")
    } else if (!"group" %in% colnames(provoc_obj)) {
        if(ncol(provoc_obj) > 4) 
            warning("Detected extra information, but plotting results as if they're a single sample.")
    
        gg <- gg +
            aes(x = date, y = rho, fill = variant) +
            labs(x = "Proportion", y = NULL, fill = "Lineage")
    } else {
        gg <- gg +
            aes(x = group, y = rho, fill = variant) +
            labs(x = "Proportion", y = NULL, fill = "Lineage")
    }
    print(gg)
}

#' Extract the variant matrix used to fit the model
#' @export
get_varmat <- function(provoc_obj) {
    attributes(provoc_obj)$variant_matrix
}


#' Extract just the results of lineage estimation
#' @export
get_res <- function(provoc_obj) {
    as.data.frame(provoc_obj)
}

#' Extract convergence information
#' @export
get_convergence <- function(provoc_obj) {
    attributes(provoc_obj)$convergence
}
