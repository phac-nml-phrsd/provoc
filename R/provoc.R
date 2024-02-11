#' Proportions of Variants of Concern (provoc) Analysis
#'
#' Applies a binomial GLM to analyze COVID-19 Variant of Concern proportions. 
#' It allows flexible lineage and mutation definitions.
#'
#' @param formula A formula for the binomial model, like cbind(count, coverage) ~ . 
#' @param data Data frame containing count, coverage, and lineage columns.
#' @param mutation_defs Optional mutation definitions; if NULL, uses astronomize().
#' @param by Column name to group and process data.
#' @param update_interval Interval for progress messages (0 to suppress).
#' @param verbose TRUE to print detailed messages.
#'
#' @return An object of class 'provoc' with GLM results per group.
#'
#' @examples
#' Example using the 'Baaijens' dataset and 'astronomize' mutation definitions.
#' library(provoc)
#'
#' # Load the Baaijens dataset
#' data("Baaijens")
#'
#' # In this example mutation_defs is NULL, so the default 'astronomize' function is used
#'
#' Baaijens$mutation <- parse_mutations(Baaijens$label)
#'
#' # Fit the model using 'provoc'
#' res <- provoc(formula = cbind(count, coverage) ~ B.1.1.7 + B.1.617.2,
#'               data = Baaijens, by = "sample_id", mutation_defs = NULL)
#'
#' # Check convergence
#' print(convergence(res))
#'
#' @export
provoc <- function(formula, data, mutation_defs = NULL, by = NULL, update_interval = 20, verbose = TRUE) {
    # Validate inputs
    if (!inherits(formula, "formula")) {
        stop("Argument 'formula' must be a formula.")
    }
    if (!is.data.frame(data)) {
        stop("Argument 'data' must be a data frame.")
    }

    # Use the astronomize function if mutation_defs is NULL
    if (is.null(mutation_defs)) {
        mutation_defs <- provoc:::astronomize()
    }

    # Prepare data based on formula and mutation_defs
    response_vars <- all.vars(formula)[1:2]  # Extract count and coverage column names from formula
    lineage_vars <- setdiff(names(mutation_defs), response_vars)

    # Check for required columns in data
    if (!all(c(response_vars, lineage_vars) %in% names(data))) {
        stop("Not all required columns found in 'data'.")
    }

    # Process the formula
    model_formula <- as.formula(paste("cbind(", paste(response_vars, collapse = ", "), ") ~ ."))

    # Conditional data grouping
    if (is.null(by)) {
        # Treat entire dataset as a single group
        grouped_data <- list(all_data = data)
    } else {
        if (!by %in% names(data)) {
            stop("Column specified in 'by' not found in data.")
        }
        grouped_data <- split(data, data[[by]])
    }

    # Initialize results
    res_list <- vector(mode = "list", length = length(grouped_data))
    names(res_list) <- names(grouped_data)

    # Group counter i that have been processed / total groups
    i <- 0

    # Process each group or the entire dataset
    for (group_name in names(grouped_data)) {
        group_data <- grouped_data[[group_name]]

        # Prepare coco df for provoc_optim
        coco <- group_data[, c("count", "coverage", "mutation")]
        varmat <- mutation_defs
        
        optim_results <- provoc:::provoc_optim(coco = coco, varmat = varmat)
        
        # Process optim_results to extract data for res_list
        point_est <- optim_results$res_df
        convergence <- optim_results$convergence

        # Extract relevant results
        point_est <- optim_results$res_df
        convergence <- optim_results$convergence
        
        # Store results
        res_list[[group_name]] <- list(point_est = point_est, convergence = convergence)
        
        # Increment counter and display progress message
        i <- i + 1
        if (verbose && update_interval > 0 && !((i - 1) %% update_interval)) {
            message(sprintf("Processed %d / %d groups", i, length(grouped_data)))
        }
    }

    # Combine results
    final_results <- do.call(rbind, res_list)

    # Ensure object is of class 'provoc'
    class(final_results) <- "provoc"

    return(final_results)
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

    if(any(!as.logical(conv$convergence))) {
        if(verbose) print(conv[which(!as.logical(conv$convergence)), -2])
        return(invisible(FALSE))

    } else {
        if(verbose) cat("All samples converged\n")
        return(invisible(TRUE))

    }
}


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
#' predicted_results <- predict.provoc(provoc_obj)
predict.provoc <- function(provoc_obj, newdata = NULL, type = NULL, se.fit = NULL,
                    dispersion = NULL, terms = NULL, na.action = NULL) {
                        
    # Ensure that provoc_obj is a list and not NULL
    if (is.null(provoc_obj) || !is.list(provoc_obj)) {
        stop("Input 'provoc_obj' must be a non-null list.")
    }
    
    # Check if the required elements are present in the input
    if (!all(c("proportions", "variant_matrix") %in% names(provoc_obj))) {
        stop("Input 'provoc_obj' must contain 'proportions' and 'variant_matrix'.")
    }

    proportions <- provoc_obj$proportions
    variant_matrix <- provoc_obj$variant_matrix

    # Check for numeric types and non-empty inputs
    if (!is.numeric(proportions) || length(proportions) == 0) {
        stop("'proportions' must be a non-empty numeric vector.")
    }
    if (!is.matrix(variant_matrix) || any(dim(variant_matrix) == 0)) {
        stop("'variant_matrix' must be a non-empty numeric matrix.")
    }

    # Error checking for dimensions
    if (ncol(variant_matrix) != length(proportions)) {
        stop("Dimensions of 'variant_matrix' and 'proportions' do not match.")
    }
    
    # Perform matrix multiplication for predictions
    results <- variant_matrix %*% proportions
    
    return(results)  # Return results in the same order
}
