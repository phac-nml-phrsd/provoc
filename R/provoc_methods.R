
#' Predict using Proportions of Variants of Concern
#'
#' Takes a named list with an estimate of the proportions and the
#' associated variant matrix, performs matrix multiplication to
#' predict outcomes, and returns results in the same order as the original data.
#'
#' @param provoc_obj Named list with `proportions` and `variant_matrix`.
#' @param newdata Not yet implemented.
#' @param type Not yet implemented.
#' @param se.fit Not yet implemented.
#' @param dispersion Not yet implemented.
#' @param terms Not yet implemented.
#' @param na.action Not yet implemented.
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

    proportions <- as.numeric(provoc_obj$rho)
    variant_matrix <- get_mutation_defs(provoc_obj)
    if (any(!rownames(variant_matrix) %in% provoc_obj$variant)) {
        stop("Variant matrix does not match variants in results")
    }

    results <- proportions %*% variant_matrix[provoc_obj$variant, ]

    return(results)
}

resids.provoc <- function(provoc_obj, type = "deviance") {
    # TODO: Calculate deviance residuals
}
