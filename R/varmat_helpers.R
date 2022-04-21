#' Create a variant matrix from a list of variants with their mutations
#' 
#' The input should be a named list where the names represent lineages and the values are vectors of mutations. The output will be a valid variant matrix.
#' 
#' @param variant_list A named list of vectors of mutations. Mutation names should be in the same format as those in your data, otherwise post-processing will be required.
#' 
#' @return A variant matrix (rownames are variants, colnames are mutations, entry i,j is 1 if variant i contains mutation j, 0 otherwise).
#' @export
varmat_from_list <- function(variant_list) {
    if (!is.list(variant_list)) {
        stop("Input must be a list")
    } else if (is.null(names(variant_list))) {
        stop("Input must be a named list.")
    } else if (any(sapply(variant_list) != "character")) {
        stop("Input must be a named list of character vectors")
    }

    # bind_rows forces all columns to be the same, adds NAs
    varmat <- dpylyr::bind_rows(lapply(variant_list, function(x) {
        # a data frame of 1s
        res <- as.data.frame(matrix(1, nrow = 1, ncol = length(x)))
        names(res) <- x
        res
    }))

    # Replace NAs induced by bind_rows
    varmat[is.na(varmat)] <- 0
    row.names(varmat) <- names(variant_list)
    varmat <- as.matrix(varmat)
    varmat
}

#' Given a vector of lineage names, uses GenBank data to determine the variant matrix.
#' 
#' Requires one of a few potentially arbitrary choices: Use the \code{max_n} most common mutations for a given variant or use the mutations in the top \code{top_quantile} quantile.
#' 
#' See \code{?mutations_by_lineage} for more information about the GenBank data used and see \code{unique(mutations_by_lineage$lineage)} to see available lineage names. 
#' 
#' @param variants A character vector of variant names in the same format as the 
#' @param max_n The maximum number of mutations from a given lineage to retain.
#' @param top_quantile Only take mutations that are in the top \code{top_quantile} quantile. E.g. 0.05 gives the mutations in more than 95% of the sequences of that lineage.
#' 
#' @return A variant matrix (rownames are variants, colnames are mutations, entry i,j is 1 if variant i contains mutation j, 0 otherwise).
#' @export
varmat_from_variants <- function(variants, 
    max_n = NULL, top_quantile = NULL) {
    if(is.null(max_n) & is.null(top_quantile)) {
        warning("Parameters not set, using all mutations (varmat will be massive and estimation will be slow and fraught with correlations)")
        top_quantile <- 0
        max_n <- 700
    } else if (is.null(top_quantile)) {
        top_quantile <- 0
    } else if (is.null(max_n)) {
        max_n <- 700
    }

    variant_list <- lapply(variants, function(x) {
        m <- mutations_by_lineage[mutations_by_lineage$lineage == x,]
        m$mutation[m$count <= max_n & m$count < quantile(m$count, 1 - top_quantile)]
    })
    names(variant_list) <- variants

    varmat_from_list(variant_list)
}


#' Use mutations found in the data to determine the variant matrix
#' 
#' Uses the included GenBank mutation list (\code{?mutations_by_lineage}) to determine a variant matrix that is guaranteed to have the mutations from your data (assuming you provide the correct mutation name format).
#' 
#' @param type A vector containing either +, -, or ~, representing insertions, deletions, and mutations, respectively. Can also be i, d, or m. 
#' @param pos A vector containing the 1-indexed positions.
#' @param alt A vector containing the alternative (nucleotide(s) for mutations and insertions or the number of deletions).
#' @param aa The vector of mutations in amino acid format (like those create by \code{parse_mutation}). TODO: Not yet implemented.
#' @param data A data frame containing columns labelled "type", "pos", and "alt", each with entries in the expected format. TODO: Not yet implemented.
#' @param mutation_format "tpa" (default) for \code{type|pos|alt} (e.g. \code{~|2832|G}; note that it includes the pipes "|") or "aa" for amino acid (e.g. \code{aa:orf1a:K856R}).
#' @param max_n The maximum number of mutations from a given lineage to retain.
#' @param top_quantile Only take mutations that are in the top \code{top_quantile} quantile. E.g. 0.95 gives the mutations that are observed in more than 95% of the sequences of that lineage.
#' 
#' @return A variant matrix (rownames are variants, colnames are mutations, entry i,j is 1 if variant i contains mutation j, 0 otherwise).
#' @export
varmat_from_data <- function(type = NULL, pos = NULL, alt = NULL, 
    aa = NULL, data = NULL, 
    mutation_format = c("tpa", "aa"),
    max_n = NULL, top_quantile = NULL) {

    # TODO: Error checking (same size, not all null, what if type and aa are both specified, etc.)
    if(is.null(max_n) & is.null(top_quantile)) {
        warning("Parameters not set, using all mutations")
        top_quantile <- 1
        max_n <- 700
    } else if (is.null(top_quantile)) {
        top_quantile <- 1
    } else if (is.null(max_n)) {
        max_n <- 700
    }

    mutations2 <- mutations_by_lineage
    if("m" %in% type) {
        mutations2$mutation <- gsub("~", "m", mutations2$mutation)
        mutations2$mutation <- gsub("\\+", "i", mutations2$mutation)
        mutations2$mutation <- gsub("-", "d", mutations2$mutation)
    }


    
    mutations <- unique(paste(type, pos, alt, sep = "|"))
    varmat <- dplyr::bind_rows(lapply(unique(mutations_by_lineage$lineage), function(x) {
        m <- mutations_by_lineage[mutations_by_lineage$lineage == x,]
        m <- m[m$count <= max_n & m$count <= quantile(m$count, top_quantile), ]
        mutations_present <- m$mutation[m$mutation %in% mutations]
        # If no mutations, create a dummy dataframe with two 
        # columns and column names that can be removed.
        if(length(mutations_present) == 0) {
            mutations_present <- c("None", "Found")
        }
        res <- as.data.frame(matrix(1, nrow = 1, ncol = length(mutations_present)))
        names(res) <- mutations_present
        res
    }))

    # TODO: Warn about unused mutations

    varmat <- varmat[, !names(varmat) %in% c("None", "Found")]
    varmat[is.na(varmat)] <- 0
    varmat <- as.matrix(varmat)
    rownames(varmat) <- unique(mutations_by_lineage$lineage)
    varmat <- varmat[rowSums(varmat) > 0, ]
    if (mutation_format[1] == "aa") {
        for(col in 1:ncol(varmat)) {
            tpa <- strsplit(colnames(varmat)[col], split = "\\|")[[1]]
            colnames(varmat)[col] <- parse_mutation(tpa[1], tpa[2], tpa[3])
        }
    }
    varmat
}
