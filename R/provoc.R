#' Proportions of Variants of Concern (provoc) Analysis
#'
#' Applies provoc_optim to analyze COVID-19 Variant of Concern proportions.
#' It allows flexible lineage and mutation definitions.
#'
#' @param formula A formula for the binomial model, like cbind(count, coverage) ~ .
#' @param data Data frame containing count, coverage, and lineage columns.
#' @param mutation_defs Optional mutation definitions; if NULL, uses astronomize().
#' @param by Column name to group and process data.
#' @param update_interval Interval for progress messages (0 to suppress).
#' @param verbose TRUE to print detailed messages.
#'
#' @return Returns an object of class 'provoc' with results from applying `provoc_optim` to the input data. The object contains the following attributes:
#'  - proportions: Estimated proportions vector for each variant of concern.
#'  - variant_matrix: Mutation definitions used for analysis, provided `mutation_defs` or default `astronomize()`.
#'
#' Outputs necessary information for subsequent analysis, including the use of the `predict.provoc()`.
#'
#' @examples
#' library(provoc)
#' # Load a dataset
#' data("Baaijens")
#' # Prepare the dataset
#' Baaijens$mutation <- parse_mutations(Baaijens$label)
#'
#' # Analyze the dataset using the default mutation definitions
#' res <- provoc(formula = cbind(count, coverage) ~ B.1.1.7 + B.1.617.2, data = Baaijens, by = "sample_id")
#'
#' # Check for analysis convergence
#' print(convergence(res))
#'
#' # Use the results for prediction
#' predicted_values <- predict.provoc(res)
#'
#' @export

provoc <- function(formula, data, mutation_defs = NULL, by = NULL,
    update_interval = 20, verbose = TRUE) {
    # Initial validation and processing
    validate_inputs(formula, data)
    mutation_defs <- as.matrix(process_mutation_defs(mutation_defs))

    # Find out which column might define the mutations
    mutation_col <- NULL
    for (col in colnames(data)) {
        if (any(colnames(mutation_defs) %in% data[, col])) {
            mutation_col <- col
        }
    }
    if (is.null(mutation_col)) {
        stop("No column contains the mutations in mutation_defs.")
    }

    # Extract components from the formula
    components <- extract_formula_components(formula, data,
        mutation_defs, mutation_col, by)
    data <- components$data
    mutation_defs <- components$mutation_defs

    # Fuse data with mutation definitions
    data <- provoc:::fuse(data, mutation_defs, verbose = verbose)

    # Group the fused data for processing
    if (!is.null(by)) {
        if (!by %in% names(data)) {
            stop("Column specified in 'by' not found in data.")
        }
        grouped_data <- split(data, data[[by]])
    } else {
        # If no grouping is specified, treat the entire dataset as one group
        grouped_data <- list(all_data = data)
    }

    # Proceed with processing each group
    res_list <- process_optim(grouped_data, mutation_defs, by)

    # Combine results and ensure object is of class 'provoc'
    final_results <- do.call(rbind, res_list)
    for (i in seq_along(res_list)) {
        if (i == 1) {
            final_results <- as.data.frame(res_list[[1]])
            if (!is.null(by)) {
                final_results$group <- unique(data[, by])[1]
            }
        } else {
            res_list[[i]]$group <- unique(data[, by])[i]
            final_results <- rbind(final_results, res_list[[i]])
        }
    }
    row.names(final_results) <- NULL

    provoc_obj <- final_results
    attr(provoc_obj, "variant_matrix") <- mutation_defs
    attr(provoc_obj, "formula") <- formula
    attr(provoc_obj, "convergence") <- attributes(res_list)$convergence
    class(provoc_obj) <- c("provoc", "data.frame")

    return(provoc_obj)
}

#' Validate Inputs for provoc
#'
#' Checks if the provided formula and data frame are valid for analysis.
#'
#' @param formula [stats]{formula}, specifying the model to be fitted.
#' @param data A data frame containing the variables in the model.
#'
#' @return None
#' @examples
#' # This function is internally used and not typically called by the user.
validate_inputs <- function(formula, data) {
    if (!inherits(formula, "formula"))
        stop("Argument 'formula' must be a formula.")

    if (!is.data.frame(data)) stop("Argument 'data' must be a data frame.")
}


#' Extract Formula Components
#'
#' Extracts and processes components from the formula provided to the provoc function.
#'
#' @param formula The formula input by the user.
#' @param data The dataframe containing the dataset.
#' @param mutation_defs A matrix containing mutation definitions.
#'
#' @return A list containing `data`, a dataframe filtered based on the formula's LHS
#' and `mutation_defs`, a matrix filtered to only include mutations on the formula's RHS
#' @examples
#' This function is internally used and not typically called by the user.
extract_formula_components <- function(formula, data, mutation_defs, mutation_col, by_col) {
    # Extract LHS and RHS of the formula
    formula_str <- deparse(formula)
    formula_parts <- strsplit(formula_str, "~")[[1]]
    lhs <- trimws(formula_parts[1])
    rhs <- trimws(formula_parts[2])

    # Split RHS by '+' and trim whitespace
    variant_names <- strsplit(rhs, "\\+")[[1]]
    variant_names <- sapply(variant_names, trimws)

    # Extract necessary data based on LHS
    response_vars <- all.vars(formula[[2]])
    necessary_data <- data[, c(mutation_col, response_vars, by_col), drop = FALSE]
    colnames(necessary_data) <- c("mutation", "count", "coverage", by_col)

    # Validate and subset mutation definitions based on RHS variants
    if (!is.null(mutation_defs) && length(variant_names) > 0) {
        missing_variants <- setdiff(variant_names, rownames(mutation_defs))
        if (length(missing_variants) > 0) {
            stop("These variants from the formula are not in mutation_defs: ",
                paste(missing_variants, collapse = ", "), ".")
        }
        necessary_mutation_defs <- mutation_defs[variant_names, , drop = FALSE]
    } else {
        necessary_mutation_defs <- mutation_defs
    }

    return(list(data = necessary_data, mutation_defs = necessary_mutation_defs))
}



#' Process Mutation Definitions
#'
#' Handles mutation definitions by using the provided matrix or generating it using \code{astronomize()}.
#'
#' @param mutation_defs A matrix of mutation definitions or NULL to use the default generated by \code{astronomize()}.
#'
#' @return A matrix of mutation definitions ready for analysis.
#' @examples
#' # This function is internally used and not typically called by the user.
process_mutation_defs <- function(mutation_defs) {
    if (is.null(mutation_defs)) {
        return(provoc:::astronomize())
    }
    if (!is.matrix(mutation_defs)) {
        mutation_defs <- tryCatch(as.matrix(mutation_defs), error = function(e) e)
        if ("error" %in% class(mutation_defs)) {
            stop("mutation_defs must be a matrix or something coercible into a matrix")}
    }
    return(as.matrix(mutation_defs))
}


#' Prepare and Fuse Data
#'
#' Prepares the data based on the grouping variable and applies the \code{fuse} function.
#'
#' @param data A data frame containing the variables in the model.
#' @param mutation_defs A matrix of mutation definitions.
#' @param by An optional string specifying the column name to group the data by.
#' @param verbose
#'
#' @return A list containing the fused data frame and the grouped data as a list (if applicable).
#' @examples
#' # This function is internally used and not typically called by the user.
prepare_and_fuse_data <- function(data, mutation_defs, by, verbose) {
    if (!is.null(by)) {
        if (!by %in% names(data))
            stop("Column specified in 'by' not found in data.")

        grouped_data <- split(data, data[[by]])
    } else {
        # Treat entire dataset as a single group
        grouped_data <- list(all_data = data)
    }
    return(grouped_data)
}


#' Process Optimization
#'
#' Processes each group or the entire dataset through \code{provoc_optim} and collects results.
#'
#' @param grouped_data A list containing data frames for each group to be processed.
#' @param mutation_defs A matrix of mutation definitions.
#'
#' @return A list of results for each group, including point estimates and convergence information.
#' @examples
#' # This function is internally used and not typically called by the user.
process_optim <- function(grouped_data, mutation_defs, by) {
    res_list <- vector("list", length = length(grouped_data))
    names(res_list) <- names(grouped_data)
    convergence_list <- vector("list", length = length(grouped_data))
    names(convergence_list) <- names(grouped_data)

    for (group_name in names(grouped_data)) {
        group_data <- grouped_data[[group_name]]
        coco <- group_data[, c("count", "coverage", "mutation", by)]
        fused <- provoc:::fuse(coco, mutation_defs, verbose = FALSE)
        fissed <- fission(fused)
        coco <- fissed$coco
        varmat <- fissed$varmat

        optim_results <- provoc:::provoc_optim(coco = coco, varmat = varmat)
        res_list[[group_name]] <- optim_results$res_df
        convergence_list[[group_name]] <- optim_results$convergence
    }
    attr(res_list, "convergence") <- convergence_list
    return(res_list)
}


#' Check if provoc converged
#'
#' If converged, returns True and prints a message. Otherwise, prints the samples and the note giving hints as to why it didn't converge.
#'
#' @param res The result of \code{provoc()}
#' @param verbose Print a message to the screen?
#'
#' @return Invisbly returns TRUE if all samples converged, false otherwise.
#' @export
convergence <- function(res, verbose = TRUE) {
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
