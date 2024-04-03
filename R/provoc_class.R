#' Print the results of lineage abundance estimation
#'
#' @inherit summary.provoc
#' @param provoc_obj The resulting object of class provoc to be used in printing.
#' @param n The number of rows of the results dataframe to print
#'
#' @export
print.provoc <- function(provoc_obj, n = 6) {
    cat("Call: ")
    print(attributes(provoc_obj)$formula)
    cat("\n\n")
    all_conv <- unlist(attributes(provoc_obj)$convergence)
    if (any(!all_conv)) {
        cat("Some models did not converge:\n")
        print(unlist(attributes(provoc_obj)$convergence)[!all_conv])
    } else {
        cat("All models converged.\n")
    }
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
#' @param provoc_obj The result of provoc(), or an object coerced via as.provoc().
#'
#' @export
summary.provoc <- function(provoc_obj) {
    cat("Call: ")
    print(attributes(provoc_obj)$formula)
    cat("\n")
    fitted <- predict(provoc_obj)

    bootstrap_cor <- attributes(provoc_obj)$bootstrap_cor

    cat(summarise_variants(provoc_obj))
    # TODO: Data summary
    # Number of mutations used in the fitting,
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
    return(invisible(list(formula = formula, resids = fitted,
                boot_cor = bootstrap_cor)))
}

#' Print the summary of a provoc object
print.summary.provoc <- function(summary.provoc) {
    # TODO: All of it.
}

#' Plot the results of model fitting
#'
#' @inherit summary.provoc
#' @param provoc_obj Resulting object of class provoc to be used in plotting.
#' @param plot_type Currently only "barplot" is implemented. Residual plots and other diagnostics are works in progress.
#'
#' @export
plot.provoc <- function(provoc_obj, plot_type = c("barplot")) {
    #plot_types <- c("b", "r")
    #if(!all(plot_type %in% plot_types)) stop("Invalid plot choice.")

    mfrow <- switch(as.character(length(plot_type)),
        "1" = c(1, 1),
        "2" = c(1, 2),
        c(2, 2))
    par(mfrow = mfrow)

    # Barplot
    if (1 %in% plot_type || any(startsWith(plot_type, "b"))) {
        barplot(
            height = matrix(provoc_obj$rho,
                ncol = ifelse("group" %in% names(provoc_obj),
                    length(unique(provoc_obj$group)),
                    1)),
            names.arg = unique(provoc_obj$group),
            col = seq_along(unique(provoc_obj$variant)),
            horiz = TRUE,
            xlim = c(0, 1),
            xlab = "Proportion",
            ylab = NULL,
            las = 1)
        legend("topright",
            legend = unique(provoc_obj$variant),
            col = seq_along(unique(provoc_obj$variant)),
            pch = 15)
    }
}

#' Plot a provoc object using ggplot2
#'
#' Plots the results of estimating wastewater prevalence of SARS-CoV-2. Optionally plots the results over time if given a date column.
#'
#' @inherit summary.provoc
#' @param provoc_obj Resulting object of class provoc to be used in plotting.
#' @param date_col Optional - if there's a date column, the results are plotted over time. This can be problematic if there are multiple samples at each time point.
#'
#' @importFrom ggplot2 autoplot
#' @export
autoplot.provoc <- function(provoc_obj, date_col = NULL) {
    if (!"ggplot2" %in% .packages()) {
        stop("Please load ggplot2 before using this function.")
    }

    gg <- ggplot(provoc_obj) +
        geom_bar(stat = "identity", position = "stack") +
        lims(y = c(0, 1))
    if (!is.null(date_col)) {
        if (!inherits(provoc_obj[, date_col], "Date"))
            stop("Supplied date column does not include Date values. \nTry lubridate::ymd().")
        gg <- gg  +
            aes(x = date, y = rho, fill = variant, group = group) +
            labs(y = "Proportion", x = "Date", fill = "Lineage") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

    } else if (!"group" %in% colnames(provoc_obj)) {
        if (ncol(provoc_obj) > 4)
            warning("Detected extra information, but plotting results as if they're a single sample.")
        gg <- gg +
            aes(x = 1, y = rho, fill = variant) +
            labs(y = "Proportion", x = NULL, fill = "Lineage") +
            scale_x_continuous(minor_breaks = NULL, breaks = NULL) +
            coord_flip()

    } else {
        gg <- gg +
            aes(x = group, y = rho, fill = variant) +
            labs(y = "Proportion", x = "Group", fill = "Lineage") +
            scale_x_discrete(breaks = sort(unique(provoc_obj$group))) +
            coord_flip()

    }
    return(gg)
}

#' Extract the mutation definitions used to fit the model
#'
#' @inherit summary.provoc
#' @param provoc_obj Object of class provoc to be used to extract the mutation definitions.
#'
#' @export
get_mutation_defs <- function(provoc_obj) {
    attributes(provoc_obj)$variant_matrix
}


#' Extract just the results of lineage estimation
#'
#' @inherit summary.provoc
#' @param provoc_obj Object of class provoc to be used to extract the results.
#'
#' @export
get_res <- function(provoc_obj) {
    as.data.frame(provoc_obj)
}


#' Check if provoc converged
#'
#' If converged, returns True and prints a message. Otherwise, prints the samples and the note giving hints as to why it didn't converge.
#'
#' @param res The result of object provoc
#' @param verbose Print a message to the screen?
#'
#' @return Invisbly returns TRUE if all samples converged, false otherwise.
#' @export
get_convergence <- function(res, verbose = TRUE) {
    if (!"convergence" %in% attributes(attributes(res))$names) {
        stop("Not a result of provoc - does not have correct attributes")
    }

    conv <- attr(res, "convergence")

    if (any(!as.logical(conv$convergence))) {
        if (verbose) print(conv[which(!as.logical(conv$convergence)), -2])
        return(invisible(FALSE))

    } else {
        if (verbose) cat("All samples converged\n")
        return(invisible(TRUE))
    }
}


#' Summarise the similarities in variant matrices
#'
#' @inherit summary.provoc
summarise_variants <- function(provoc_obj) {
    similarities <- attributes(provoc_obj)$internal_data |>
        provoc:::variants_similarity() |>
        provoc:::simplify_similarity()

    msg <- ""
    if (length(similarities$Differ_by_one_or_less) > 0) {
        msg <- paste0(msg,
            "At least one pair of variants has a single difference. ",
            collapse = " ")
    }
    if (length(similarities$Jaccard_similarity) > 0) {
        msg <- paste0(msg,
            "At least one pair of variants has a Jaccard similarity > 0.99",
            collapse = " ")
    }
    if (length(similarities$Is_subset) > 0) {
        msg <- paste0(msg,
            "At least one variant is a subset of another.",
            collapse = " ")
    }
    if (length(similarities$Is_almost_subset) > 0) {
        msg <- paste0(msg,
            "At least one variant is almost a subset of another.",
            collapse = " ")
    }

    if (nchar(msg) > 0) {
        msg <- paste0(msg, "See variants_similarity() for more info.\n",
            collapse = "\n")
    } else {
        msg <- "No issues detected for mutation definition.\n"
    }

    msg
}

#' Plot the residuals, by variant
#'
#' @param provoc_obj Result of fitting provoc().
#' @param type Deviance or raw residuals.
#'
plot_resids <- function(provoc_obj, type = "deviance", by_variant = TRUE) {
    data <- attributes(provoc_obj)$internal_data
    vardf <- data[, startsWith(colnames(data), "var_")]
    varnames <- colnames(vardf)[startsWith(colnames(vardf), "var_")]
    vardf$fitted <- as.numeric(predict(provoc_obj))
    dcov <- ifelse(data$coverage == 0, 1, data$coverage)
    vardf$residuals <- provoc:::resids.provoc(provoc_obj, type = type)

    plot(NA,
        xlim = c(0, 1),
        ylim = range(vardf$residuals, na.rm = TRUE),
        xlab = "Fitted",
        ylab = paste0(tools::toTitleCase(type), " Residuals"),
        main = "Residuals versus Fitted"
    )
    abline(h = 0, col = "lightgrey", lty = 2, lwd = 2)
    if (by_variant) {
        for (i in seq_len(ncol(vardf) - 2)) {
            points(residuals ~ fitted, data = vardf[vardf[, i] == 1, ], col = i, pch = 16)
        }
        legend("bottomright",
            legend = gsub("var_", "", varnames),
            col = seq_along(varnames),
            pch = 1)
    } else {
        points(residuals ~ fitted, data = vardf, pch = 16)
    }
}