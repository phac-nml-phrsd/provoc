
#' Predict using Proportions of Lineages
#'
#' Takes a named list with an estimate of the proportions and the
#' associated lineage matrix, performs matrix multiplication to
#' predict outcomes, and returns results in the same order as the original data.
#'
#' @param provoc_obj Named list with `proportions` and `lineage_matrix`.
#' @param newdata Not yet implemented.
#' @param type Not yet implemented.
#' @param dispersion Not yet implemented.
#' @param terms Not yet implemented.
#' @return Predicted values in the same order as the input data.
#' @export
#' @examples
#' predicted_results <- predict(provoc_obj)
predict.provoc <- function(provoc_obj,
    newdata = NULL, type = NULL,
    dispersion = NULL, terms = NULL) {

    if (!"provoc" %in% class(provoc_obj)) {
        stop("Object must be of class 'provoc'")
    }

    # Process internal data
    internal_data <- attributes(provoc_obj)$internal_data
    by_col <- attributes(provoc_obj)$by_col
    data_groups <- internal_data[, by_col]
    if (!is.null(dim(data_groups))) {
        data_groups <- rep(1, nrow(internal_data))
    }

    # Get matching lineage names
    lin_index <- startsWith(names(internal_data), "lin_")
    lin_names <- names(internal_data)[lin_index]

    # Reshape res and prepare for element-wise multiplication
    res_wide <- tidyr::pivot_wider(
        provoc_obj[, c("rho", "lineage", "group")],
        values_from = rho, names_from = lineage,
        names_prefix = "lin_"
    )
    res_wide <- res_wide[match(data_groups, res_wide$group), ]

    # Multiply the correct rows together, return as result
    rowSums(as.matrix(res_wide[, lin_names]) *
            as.matrix(internal_data[, lin_names]))

}

#' Calculate the residuals of a provoc object
#' 
#' @param provoc_obj The result of provoc()
#' @param type "deviance" or "raw"
#' 
#' @export
resid.provoc <- function(provoc_obj, type = "deviance") {
    counts <- attributes(provoc_obj)$internal_data$count
    covs <- attributes(provoc_obj)$internal_data$coverage
    covs <- ifelse(covs == 0, yes = 1, no = covs)
    preds <- predict(provoc_obj)

    raw_resids <- counts / covs - preds

    if (startsWith(x = type, prefix = "d")) {
        saturated <- dbinom(counts, size = covs,
            prob = counts / covs, log = TRUE)
        modelled <- dbinom(counts, size = covs, prob = preds)
        return(sign(raw_resids) * 2 * (saturated - modelled))
    } else {
        return(raw_resids)
    }
}

#' Calculate the residuals of a provoc object
#' 
#' @param provoc_obj The result of provoc()
#' @param type "deviance" or "raw"
#' 
#' @export
residuals.provoc <- resid.provoc
