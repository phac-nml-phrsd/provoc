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
    cat("\n")
    minfo <- attributes(provoc_obj)$mutation_info
    cat("Mutations in lineage definitions: ",
        ncol(attributes(provoc_obj)$variant_matrix),
        "\n")
    if (length(minfo[[1]]) <= 10) {
        cat("Mutations used in analysis/mutations in data:\n")
        cat(paste(minfo[[1]], minfo[[2]], sep = "/", collapse = "\t"))
        cat("\n\n")
    } else {
        cat("Summary of percent of mutations in data used:\n")
        print(summary(minfo[[1]] / minfo[[2]]))
        cat("\n")
    }
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
    cat("\nCall:\n")
    print(attributes(provoc_obj)$formula)
    cat("\n")

    # Deviance is intentionally misspelled by Devan
    devance_res <- residuals(provoc_obj, type = "deviance")
    cat("Deviance Residuals:\n")
    print(summary(devance_res))

    minfo <- attributes(provoc_obj)$mutation_info
    cat("\nMutations in lineage definitions:",
        ncol(attributes(provoc_obj)$variant_matrix),
        "\n")
    if (length(minfo[[1]]) <= 10) {
        cat("Mutations used in analysis/mutations in data:\n")
        cat(paste(minfo[[1]], minfo[[2]], sep = "/", collapse = "\t"))
        cat("\n")
    } else {
        cat("Summary of percent of mutations in data used:\n")
        print(summary(minfo[[1]] / minfo[[2]]))
        cat("\n")
    }

    cat("\nCoefficients:\n")
    coef_table <- as.data.frame(provoc_obj)
    coef_table_length <- seq_len(min(30,
            nrow(coef_table)))
    by_col <- attributes(provoc_obj)$by_col
    if (is.null(by_col)) {
        coef_table <- coef_table[, 1:4]
    } else {
        coef_table <- coef_table[, 1:5]
    }
    print(coef_table[
            coef_table_length, ])

    cat("\nCorrelation of coefficients:\n")
    bootstraps <- attributes(provoc_obj)$bootstrap
    if (length(bootstraps) > 1) {
        cat("\nMultiple samples detected, see boot_corr() for info.\n")
    } else {
        boot_corr <- cor(bootstraps[[1]])
        if (ncol(boot_corr) <= 6) {
            print(boot_corr)
        } else {
            cat("Top 6 are shown; see coef_cor() for more.")
            boot_corr[lower.tri(boot_corr)] <- 0
            diag(boot_corr) <- 0
            max_corrs <- apply(boot_corr, 2, max)
            top_6 <- names(sort(-max_corrs))[1:6]
            print(boot_corr[top_6, top_6])
        }
    }
    return(
        invisible(
            list(
                formula = formula,
                resids = devance_res,
                coef = coef_table,
                boot_corr = boot_corr
            )
        )
    )
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
    similarities <- attributes(provoc_obj)$similarities |>
        provoc:::simplify_similarity()

    msg <- ""
    if (length(similarities$Differ_by_one_or_less) > 0) {
        msg <- paste0(msg,
            "At least one pair of variants has a single difference.\n",
            collapse = " ")
    }
    if (length(similarities$Jaccard_similarity) > 0) {
        msg <- paste0(msg,
            "At least one pair of variants has a Jaccard similarity > 0.99.\n",
            collapse = " ")
    }
    if (length(similarities$is_subset) > 0) {
        msg <- paste0(msg,
            "At least one variant is a subset of another.\n",
            collapse = " ")
    }
    if (length(similarities$is_almost_subset) > 0) {
        msg <- paste0(msg,
            "At least one variant is almost a subset of another.\n",
            collapse = " ")
    }

    if (nchar(msg) > 0) {
        msg <- paste0(msg, "See variants_similarity(res) for more info.\n",
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
#' @export
plot_resids <- function(provoc_obj, type = "deviance", by_variant = TRUE) {
    data <- attributes(provoc_obj)$internal_data
    vardf <- data[, startsWith(colnames(data), "var_")]
    varnames <- colnames(vardf)[startsWith(colnames(vardf), "var_")]
    vardf$fitted <- as.numeric(predict(provoc_obj))
    vardf$residuals <- provoc:::resid.provoc(provoc_obj, type = type)

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
            points(residuals ~ fitted,
                data = vardf[vardf[, i] == 1, ],
                col = i, pch = 16)
        }
        legend("bottomright",
            legend = gsub("var_", "", varnames),
            col = seq_along(varnames),
            pch = 16)
    } else {
        points(residuals ~ fitted, data = vardf, pch = 16)
    }
}

#' plot the similarities of variants
#' @export
plot_variants <- function(provoc_obj,
    type = "Jaccard_similarity", labels = TRUE) {
    
    old_par <- par()
    similarities <- attributes(provoc_obj)$similarities[[type]]
    similarities[upper.tri(similarities)] <- NA
    diag(similarities) <- NA
    similarities <- similarities[-1, -ncol(similarities)]
    rownames <- gsub("var_", "", rownames(similarities))
    colnames <- gsub("var_", "", colnames(similarities))
    atseq <- (seq_along(rownames) - 1) / (length(rownames) - 1)

    omar <- par()$mar
    par(mar = c(6, 2, 2, 6))
    image(similarities, xaxt = "n", yaxt = "n", bty = "n",
        zlim = c(0, 1), main = type,
        col = hcl.colors(
            n = length(atseq)^4,
            palette = "Spectral",
            rev = TRUE))
    axis(side = 1, at = atseq, las = 2,
        labels = rownames)
    axis(side = 4, at = atseq,
        labels = colnames, las = 1)

    if (labels) {
        sims <- round(as.numeric(similarities), 3)
        sims[is.na(sims)] <- ""
        atseqy <- rep(atseq, each = length(atseq))
        atseqx <- rep(atseq, times = length(atseq))
        text(x = atseqx, y = atseqy, labels = sims)
    }
    par(mar = omar)
}
