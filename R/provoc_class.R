#' Print the results of lineage abundance estimation
#' @export
print.provoc <- function(provoc_obj, n = 6) {
    cat("Convergence: ")
    print(!is.null(attributes(provoc_obj)$convergence))
    cat("\n")
    cat("Top", n, "variants:\n")
    print.data.frame(provoc_obj[order(-provoc_obj$rho), ][1:n, ])
}

summary.provoc <- function(provoc_obj) {
    cat("Call: ")
    print(attributes(provoc_obj)$formula)
    cat("\n")
    fitted <- predict(provoc_obj)

    return(list(formula = formula, resids = resids))
}

#' Extract the variant matrix used to fit the model
#' @export
get_varmat <- function(provoc_obj) {
    attributes(provoc_obj)$varmat
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
