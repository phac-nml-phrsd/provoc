#' Proportions of Variants of Concern (provoc) Analysis
#'
#' Applies provoc_optim to analyze COVID-19 lineage proportions.
#' It allows flexible lineage and mutation definitions.
#'
#' @param formula A formula for the binomial model, like cbind(count, coverage) ~ .
#' @param data Data frame containing count, coverage, and lineage columns.
#' @param lineage_defs Optional mutation definitions; if NULL, uses astronomize().
#' @param by Column name to group and process data. If included, the results will contain a column labelled "group".
#' @param bootstrap_samples The number of bootstrap samples to use.
#' @param update_interval Interval for progress messages (0 to suppress).
#' @param verbose TRUE to print detailed messages.
#' @param annihilate TRUE to remove duplicate lineages from the data
#'
#' @return Returns an object of class 'provoc' with results from applying `provoc_optim` to the input data. The object has methods for print(), plot(), and summary(), but otherwise is treated as a dataframe. Fit information can be extracted as follows:
#'  - convergence(res)
#'  - get_lineage_defs(res): Mutation definitions used for analysis, provided `lineage_defs` or default `astronomize()`.
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
#' print(get_convergence(res))
#'
#' # Use the results for prediction
#' predicted_values <- predict.provoc(res)
#'
#' @export
provoc <- function(formula, data, lineage_defs = NULL, by = NULL,
    bootstrap_samples = 0, update_interval = 20,
    verbose = FALSE, annihilate = FALSE) {
    #creating original copy of data for later use
    data_copy <- data

    # Initial validation and processing
    validate_inputs(formula, data)
    lineage_defs <- as.matrix(process_lineage_defs(lineage_defs))

    # Find out which column might define the mutations
    find_muation_column_output <- find_mutation_column(data, lineage_defs)
    lineage_defs <- find_muation_column_output[[1]]
    mutation_col <- find_muation_column_output[[2]]

    # Find how many mutations they have in common
    if (is.null(by)) {
        mutations_used <- length(intersect(
            data[, mutation_col],
            colnames(lineage_defs)))
        mutations_in_data <- length(unique(data[, mutation_col]))
    } else {
        mutations_used <- sapply(unique(data[, by]), function(x) {
            length(intersect(
                data[data[, by] == x, mutation_col],
                colnames(lineage_defs)
            ))
        })
        mutations_in_data <- sapply(unique(data[, by]),
            function(x) {
                length(unique(data[data[, by] == x, mutation_col]))
            }
        )
    }
    mutation_info <- list(mutations_used = mutations_used,
        mutations_in_data = mutations_in_data)

    # Extract components from the formula
    components <- extract_formula_components(formula, data,
        lineage_defs, mutation_col, by)
    data <- components$data
    lineage_defs <- components$lineage_defs

    # Fuse data with mutation definitions
    data <- provoc:::fuse(data, lineage_defs, verbose = verbose)

    #Finding summary statistics while lineage_defs is still in lineage_defs form
    similarities <- provoc:::lineage_similarity(data)

    #remove identical lineages
    data <- remove_identical_lineages(data, annihilate)

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
    res <- process_optim(grouped_data, lineage_defs,
        by, bootstrap_samples, verbose)
    res_list <- res$res_list

    # Combine results and ensure object is of class 'provoc'
    final_results <- do.call(rbind, res_list)
    for (i in seq_along(res_list)) {
        if (i == 1) {
            final_results <- as.data.frame(res_list[[1]])
            if (!is.null(by)) {
                final_results$group <- unique(data[, by])[1]
            } else {
                final_results$group <- 1
            }
        } else {
            res_list[[i]]$group <- unique(data[, by])[i]
            final_results <- rbind(final_results, res_list[[i]])
        }
    }
    row.names(final_results) <- NULL

    # Find columns that are constant to by and add them to necessary data
    constant_columns <- constant_with_by(data_copy, by)
    if (!is.null(constant_columns)) {
        final_results <- dplyr::left_join(
            final_results,
            constant_columns,
            by = c("group" = "sra"))
    }

    provoc_obj <- final_results
    attr(provoc_obj, "lineage_defs") <- lineage_defs
    attr(provoc_obj, "formula") <- formula
    attr(provoc_obj, "convergence") <- res$convergence_list
    attr(provoc_obj, "bootstrap") <- res$boot_list
    attr(provoc_obj, "similarities") <- similarities
    attr(provoc_obj, "internal_data") <- data
    attr(provoc_obj, "by_col") <- by
    attr(provoc_obj, "mutation_info") <- mutation_info

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

#' Finds the right mutation column in lineage_defs.
#'
#' Checks the columns of lineage_defs and if not found in columns will check the rows
#' and transpose lineage_defs.
#'
#' @param data Data frame containing count, coverage, and lineage columns
#' @param lineage_defs Optional mutation definitions
#'
#' @return A list of length 2 containing the possibly transposed lineage_defs and
#' the right mutation column.
find_mutation_column <- function(data, lineage_defs) {
    mutation_col <- NULL
    for (col in colnames(data)) {
        if (any(colnames(lineage_defs) %in% data[, col])) {
            mutation_col <- col
            return(list(lineage_defs, mutation_col))
        } else if (any(rownames(lineage_defs) %in% data[, col])) {
            mutation_col <- col
            return(list(t(lineage_defs), mutation_col))
        }
    }
    if (is.null(mutation_col)) {
        stop("No column contains the mutations in lineage_defs.")
    }
}

#' Remove Identical Lineages
#'
#' To ensure the predictor matrix is singular, if annihilate is TRUE the function returns all the mutations with a unique combination of lineages. If FALSE it warns user of duplicate lineages
#'
#' @param fused_df The fused data frame from the \code{fuse()} function
#' @param annihilate if TRUE will remove all duplicate lineages, if FALSE will warn user if there is duplicate lineages
#'
#' @return A data frame with no mutations that have a duplicate combination of lineages
#' @examples
#' # This function is internally used and not typically called by the user.
remove_identical_lineages <- function(fused_df, annihilate) {
    subset_of_lineages <- dplyr::select(fused_df, contains("lin_"))
    if (annihilate) {
        unique_subset_of_lineages <- t(unique(t(subset_of_lineages)))
        if (ncol(subset_of_lineages) != ncol(unique_subset_of_lineages)) {
            print(paste0("Lineage ",
                colnames(subset_of_lineages)[
                    !(colnames(subset_of_lineages) %in%
                            colnames(unique_subset_of_lineages))],
                       " is a duplicate lineage and has been removed from the dateframe"))
        }
        return(cbind(dplyr::select(fused_df,
            !contains("lin_")), unique_subset_of_lineages))
    } else {
        names_of_duplicate_lin <- names(which(duplicated(t(fused_df))))
        for (lin in names_of_duplicate_lin) {
            other_lineages <- dplyr::select(subset_of_lineages,!contains(lin))
            duplicated_with_lin <- c(lin)
            for (lin_i in colnames(other_lineages)) {
                if (all(other_lineages[[lin_i]] == subset_of_lineages[[lin]])) {
                    duplicated_with_lin <- append(duplicated_with_lin, lin_i)
                }
            }
            warning("Lineages ",
                paste(duplicated_with_lin, collapse = ", "),
                " are duplicates of eachother")
        }
        return(fused_df)
    }
}

#' Extract Formula Components
#'
#' Extracts and processes components from the formula provided to the provoc function.
#'
#' @param formula The formula input by the user.
#' @param data The dataframe containing the dataset.
#' @param lineage_defs A matrix containing mutation definitions.
#' @param mutation_col The column containing the mutations.
#' @param by_col The column to group the data by.
#'
#' @return A list containing `data`, a dataframe filtered based on the formula's LHS
#' and `lineage_defs`, a matrix filtered to only include mutations on the formula's RHS
#' @examples
#' This function is internally used and not typically called by the user.
extract_formula_components <- function(formula, data,
    lineage_defs, mutation_col, by_col) {
    # Extract LHS and RHS of the formula
    formula_str <- deparse(formula)
    formula_parts <- strsplit(formula_str, "~")[[1]]
    lhs <- trimws(formula_parts[1])
    rhs <- trimws(formula_parts[2])

    # Split RHS by '+' and trim whitespace
    lineage_names <- strsplit(rhs, "\\+")[[1]]
    lineage_names <- sapply(lineage_names, trimws)
    if (length(lineage_names) == 1) {
        if (lineage_names == ".") {
            lineage_names <- rownames(lineage_defs)
        }
    }

    # Extract necessary data based on LHS
    response_vars <- all.vars(formula[[2]])
    necessary_data <- data[, c(mutation_col, response_vars, by_col),
        drop = FALSE]
    colnames(necessary_data) <- c("mutation", "count", "coverage", by_col)

    # Validate and subset mutation definitions based on RHS lineages
    if (!is.null(lineage_defs) && length(lineage_names) > 0) {
        missing_lineages <- setdiff(lineage_names, rownames(lineage_defs))
        if (length(missing_lineages) > 0) {
            stop("These lineages from the formula are not in lineage_defs: ",
                paste(missing_lineages, collapse = ", "), ".")
        }
        necessary_lineage_defs <- lineage_defs[lineage_names, , drop = FALSE]
    } else {
        necessary_lineage_defs <- lineage_defs
    }

    return(list(data = necessary_data, lineage_defs = necessary_lineage_defs))
}



#' Process Mutation Definitions
#'
#' Handles mutation definitions by using the provided matrix or generating it using \code{astronomize()}.
#'
#' @param lineage_defs A matrix of mutation definitions or NULL to use the default generated by \code{astronomize()}.
#'
#' @return A matrix of mutation definitions ready for analysis.
#' @examples
#' # This function is internally used and not typically called by the user.
process_lineage_defs <- function(lineage_defs) {
    if (is.null(lineage_defs)) {
        return(provoc:::astronomize())
    }
    if (!is.matrix(lineage_defs)) {
        lineage_defs <- tryCatch(as.matrix(lineage_defs),
            error = function(e) e)
        if ("error" %in% class(lineage_defs)) {
            stop("lineage_defs must be a matrix or something coercible into a matrix")}
    }
    return(as.matrix(lineage_defs))
}


#' Prepare and Fuse Data
#'
#' Prepares the data based on the grouping variable and applies the \code{fuse} function.
#'
#' @param data A data frame containing the variables in the model.
#' @param lineage_defs A matrix of mutation definitions.
#' @param by An optional string specifying the column name to group the data by.
#' @param verbose
#'
#' @return A list containing the fused data frame and the grouped data as a list (if applicable).
#' @examples
#' # This function is internally used and not typically called by the user.
prepare_and_fuse_data <- function(data, lineage_defs, by, verbose) {
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
#' @param lineage_defs A matrix of mutation definitions.
#' @param by An optional string specifying the column name to group the data by.
#' @param bootstrap_samples The number of bootstrap samples to use.
#' @param verbose TRUE to print detailed messages.
#'
#' @return A list of results for each group, including point estimates and convergence information.
#' @examples
#' # This function is internally used and not typically called by the user.
process_optim <- function(grouped_data, lineage_defs, by, bootstrap_samples, verbose) {
    res_list <- vector("list", length = length(grouped_data))
    names(res_list) <- names(grouped_data)
    convergence_list <- vector("list", length = length(grouped_data))
    names(convergence_list) <- names(grouped_data)
    boot_list <- vector("list", length = length(grouped_data))
    names(boot_list) <- names(grouped_data)

    for (group_name in names(grouped_data)) {
        group_data <- grouped_data[[group_name]]
        coco <- group_data[, c("count", "coverage", "mutation", by)]
        fused <- provoc:::fuse(coco, lineage_defs, verbose = FALSE)
        fissed <- fission(fused)
        coco <- fissed$coco
        lineage_defs <- fissed$lineage_defs

        optim_results <- provoc:::provoc_optim(
            coco = coco, lineage_defs = lineage_defs,
            bootstrap_samples = bootstrap_samples,
            verbose = verbose)
        res_list[[group_name]] <- optim_results$res_df
        convergence_list[[group_name]] <- optim_results$convergence
        boot_list[[group_name]] <- optim_results$bootstrap_samples
    }

    return(list(res_list = res_list,
            convergence_list = convergence_list,
            boot_list = boot_list))
}

#' Finds all columns of the data that are constant with the specified by group
#'
#' @param data Data frame containing count, coverage, and lineage columns.
#' @param by Column name to group and process data.
#'
#' @return A dataframe with the columns by and all that are constant with by
constant_with_by <- function(data, by) {
    if (is.null(by)) {
        return(NULL)
    }
    grouped_data <- split(data, data[[by]])
    for (group_name in names(grouped_data)) {
        group_data <- grouped_data[[group_name]]
        is_constant <- apply(group_data, 2,
            function(a) length(unique(a)) == 1)
        data <- data[, names(is_constant)[is_constant]]
    }

    return(data[!duplicated(data), ])
}
