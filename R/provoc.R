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
#' @param annihilate TRUE to remove duplicate variants from the data
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
#' res <- provoc(formula = cbind(count, coverage) ~ B.1.1.7 + B.1.617.2,
#'               data = Baaijens, by = "sample_id")
#'
#' # Check for analysis convergence
#' print(convergence(res))
#'
#' # Use the results for prediction
#' predicted_values <- predict.provoc(res)
#'
#' @export

provoc <- function(formula, data, mutation_defs = NULL, by = NULL, update_interval = 20, verbose = TRUE, annihilate = FALSE) {
    # Initial validation and processing
    validate_inputs(formula, data)
    mutation_defs <- as.matrix(process_mutation_defs(mutation_defs))
    
    # Fuse data with mutation definitions
    data <- provoc:::fuse(data, mutation_defs, verbose = verbose)
    
    #remove identical variants
    data <- remove_identical_variants(data, annihilate)
    
    # Group the fused data for processing
    if (!is.null(by)) {
        if (!by %in% names(data)) {
            stop("Column specified in 'by' not found in data.")
        }
        grouped_data <- split(data, data[[by]])
    } else {
        # If no grouping is specified, treat the entire fused dataset as a single group
        grouped_data <- list(all_data = data)
    }

    # Proceed with processing each group
    res_list <- process_optim(grouped_data, mutation_defs)
    
    # Combine results and ensure object is of class 'provoc'
    final_results <- do.call(rbind, res_list)
    class(final_results) <- "provoc"

    # 'provoc' object attributes with attributes nessessary for predict.provoc()
    return(list(
        proportions = final_results,
        variant_matrix = mutation_defs
    ))
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
    if (!inherits(formula, "formula")) stop("Argument 'formula' must be a formula.")
    if (!is.data.frame(data)) stop("Argument 'data' must be a data frame.")
}

#' Remove Identical variants
#'
#' To ensure the predictor matrix is singular, if annihilate is TRUE the function returns all the mutations with a unique combination of variants. If FALSE it warns user of duplicate variants
#'
#' @param fused_df The fused data frame from the \code{fuse()} function
#' @param annihilate if TRUE will remove all duplicate variants, if FALSE will warn user if there is duplicate variants
#'
#' @return A data frame with no mutations that have a duplicate combination of variants
#' @examples
#' # This function is internally used and not typically called by the user.
remove_identical_variants <- function(fused_df, annihilate){
    subset_of_variants <- dplyr::select(fused_df,contains("var_"))
    if (annihilate) {
        unique_subset_of_variants <- t(unique(t(subset_of_variants)))
        return(cbind(dplyr::select(fused_df, !contains("var_")), unique_subset_of_variants))
    }
    else {
        names_of_duplicate_var <- names(which(duplicated(t(fused_df))))
        for (var in names_of_duplicate_var) {
            other_variants <- dplyr::select(subset_of_variants,!contains(var))
            duplicated_with_var <- c(var)
            for (var_i in colnames(other_variants)) {
                if (all(other_variants[[var_i]] == subset_of_variants[[var]])) {
                    duplicated_with_var <- append(duplicated_with_var, var_i)
                }
            }
            warning("Variants ", paste(duplicated_with_var, collapse = ", "), " are duplicates of eachother")
        }
        return(fused_df)
    }
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
        stop("mutation_defs must be a matrix with appropriate dimension names.")
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
        if (!by %in% names(data)) stop("Column specified in 'by' not found in data.")
        grouped_data <- split(data, data[[by]])
    } else {
        grouped_data <- list(all_data = data)  # Treat entire dataset as a single group
    }
    return(list(fused_data = provoc:::fuse(data, mutation_defs, verbose = verbose), grouped_data = grouped_data))
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
process_optim <- function(grouped_data, mutation_defs) {
    res_list <- vector("list", length = length(grouped_data))
    names(res_list) <- names(grouped_data)

    for (group_name in names(grouped_data)) {
        group_data <- grouped_data[[group_name]]
        coco <- group_data[, c("count", "coverage", "mutation")]
        optim_results <- provoc:::provoc_optim(coco = coco, varmat = mutation_defs)
        res_list[[group_name]] <- list(point_est = optim_results$res_df, convergence = optim_results$convergence)
    }
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
