#' Ensure mutations are present in both coco and varmat
#'
#' Finds the intersection between the mutations present in coco and varmat. Will squash lineages together if the resulting mutation list is too similar (see details).
#'
#' @param coco A data frame with a column labelled \code{mutation}.
#' @param varmat Rownames are variants, column names are Mutations.
#' @param min_perc A variant must have at least \code{min_perc} of the mutations.
#' @param vebose Print information about mutations that were removed by the fusion. 0 (FALSE) returns errors, 1 (TRUE) returns warnings and some info about relative mutation counts, and 2 returns all mutations in each.
#'
#' @return A data frame with the same columns as coco (possibly fewer rows) and the same columns plus new columns for the variants of concern. The provoc function expects this structure.
#' @export
#'
#' @details First, the intersection of the mutations is found.
#'
#' The columns of varmat are subsetted according to this intersection. If this removes all mutations for any lineage, that lineage is removed from the study (a warning is given if the lineage was a pre-specified voc).
#'
#' After removing mutations, it's possible that some rows lose their distinctive mutations and become identical. In this case the names of the lineages are pasted together and only one of the rows are kept.
#'
#' Duplicate mutation names in coco are NOT removed. It is safe to use this function on a data frame that contains multiple samples.
fuse <- function(coco, varmat, min_perc = 0.01, verbose = FALSE) {
    if (any(colnames(coco) %in% paste0("var_", rownames(varmat)))) {
        stop("coco should not contain column names that are names of lineages. Is this object already fused?")
    }

    pre <- nrow(coco)
    shared <- intersect(unique(coco$mutation), colnames(varmat))
    # We can't say anything about mutations not in varmat.
    coco <- coco[!is.na(coco$mutation), ]
    coco <- coco[coco$mutation %in% shared, ]

    if (length(shared) < 3) {
        stop("Too few shared mutations. Are they in the same format?")
    } else if (length(shared) <= 10 && ncol(varmat) > 10 && verbose) {
        warning("Fewer than 10 shared mutations. Results may be very difficult to interpret.")
    } else if (length(shared) / pre < 0.1 && verbose) {
        warning(paste0("Less than ", length(shared) / pre,
            "% of coco's mutations are being used. Consider a larger variant matrix."))
    }
    if (verbose > 0) {
        perc_rm <- 1 - round(length(shared) / pre, 3)
        print(paste0(100 * perc_rm,
                "% of the rows of coco have been removed."))
        coco_only <- coco$mutation[!coco$mutation %in% shared]
        varmat_only <- colnames(varmat)[!colnames(varmat) %in% shared]
        print("coco-only mutations removed:")
        print(length(coco_only))
        print("varmat-only mutations removed")
        print(length(varmat_only))
    }


    # Lineages with too few mutations (less than 10% are 1s)
    # TODO: Make this a parameter?
    vari2 <- varmat[, shared]
    too_many_zeros <- apply(vari2, 1, sum) <= (ncol(vari2) * min_perc)
    vari2 <- vari2[!too_many_zeros, ]

    # Squash (nearly) identical lineages
    i <- 0
    while (i < nrow(vari2)) {
        i <- i + 1
        this_row <- vari2[i, ]
        dupes <- apply(vari2, 1, function(x) mean(this_row == x))
        dupes[i] <- 0
        if (any(dupes > 0.99)) {
            squashed <- rownames(vari2)[c(i, which(dupes > 0.99))]
            rownames(vari2)[i] <- paste(squashed, collapse = "|")
            vari2 <- vari2[dupes <= 0.99, ]
        }
    }

    vardf <- as.data.frame(t(vari2))
    names(vardf) <- paste0("var_", names(vardf))
    vardf$mutation <- rownames(vardf)

    dplyr::left_join(coco, vardf, by = "mutation")
}

#' Un-fuse coco and varmat.
#'
#' Fusion ensures that the mutation lists match and are in the correct order, but the two must be separated.
#'
#' @param fused The result of \code{fuse(coco, varmat)}
#' @param sample The name of the sample being used.
#'
#' @return A list containing coco and varmat.
fission <- function(fused, sample = NULL) {
    if (!is.null(sample)) fused <- fused[fused$sample == sample, ]
    variants <- startsWith(names(fused), "var_")
    varnames <- names(fused)[variants]
    coco <- fused[, !variants]

    vardf <- fused[, variants]
    varmat <- t(as.matrix(vardf))
    varmat <- matrix(as.numeric(varmat), ncol = ncol(varmat))
    rownames(varmat) <- gsub("var_", "", varnames)
    colnames(varmat) <- coco$mutation

    return(list(coco = coco, varmat = varmat))
}

#' Finds and prints all similarities among variants
#'
#' @param data A variant matrix
#'
#' @return A list of length 4 containing information on which variants differ by one,
#' the Jaccard similarity between variants, which variants are subsets and almost subsets
#' of each other. in is_subset and is_almost_subset a value is true if the variant of the
#' column name is a subset/almost a subset of the variant of the row name.
variants_similarity <- function(data, simplify = FALSE, almost = 1) {

    subset_of_variants <- data[, startsWith(names(data), "var_")]

    similarities <- list()

    # DIFFER BY ONE OR LESS -------------------------------
    similarities$Differ_by_one_or_less <- outer(
        colnames(subset_of_variants),
        colnames(subset_of_variants),
        function(x, y) {
            mapply(FUN = differ_by_one_or_less,
                v1 = subset_of_variants[, x],
                v2 = subset_of_variants[, y])
        })
    colnames(similarities$Differ_by_one_or_less) <- colnames(subset_of_variants)
    rownames(similarities$Differ_by_one_or_less) <- colnames(subset_of_variants)

    # JACCARD ---------------------------------------------
    similarities$Jaccard_similarity <- outer(
        colnames(subset_of_variants),
        colnames(subset_of_variants),
        function(x, y) {
            mapply(
                FUN = jaccard_simularity,
                v1 = subset_of_variants[, x],
                v2 = subset_of_variants[, y])
        })
    colnames(similarities$Jaccard_similarity) <- colnames(subset_of_variants)
    rownames(similarities$Jaccard_similarity) <- colnames(subset_of_variants)

    # SUBSET ----------------------------------------------
    similarities$is_subset <- outer(
        colnames(subset_of_variants),
        colnames(subset_of_variants),
        function(x, y) {
            mapply(
                FUN = is_subset,
                v1 = subset_of_variants[, x],
                v2 = subset_of_variants[, y])
        })
    colnames(similarities$is_subset) <- colnames(subset_of_variants)
    rownames(similarities$is_subset) <- colnames(subset_of_variants)


    # ALMOST SUBSET ---------------------------------------
    similarities$is_almost_subset <- outer(
        colnames(subset_of_variants),
        colnames(subset_of_variants),
        function(x, y) {
            mapply(
                FUN = is_almost_subset,
                v1 = subset_of_variants[, x],
                v2 = subset_of_variants[, y])
        })
    colnames(similarities$is_almost_subset) <- colnames(subset_of_variants)
    rownames(similarities$is_almost_subset) <- colnames(subset_of_variants)

    if (simplify) {
        similarities <- simplify_similarity(similarities,
            almost = almost)
    }
    return(similarities)
}

#' Simplifies variant similarity matrices for easier interpretation
#' 
#' @param similarities Result of \code{variants_similarity()}
#' @param almost Degree of similarity.
simplify_similarity <- function(similarities, almost) {


    keep_vars <- apply(X = similarities$Differ_by_one_or_less,
        MARGIN = 1,
        FUN = function(x) {
            sum(x) > 1
        })

    similarities$Differ_by_one_or_less <-
        similarities$Differ_by_one_or_less[keep_vars, keep_vars]

    diag(similarities$Jaccard_similarity) <- 0
    keep_vars <- apply(X = similarities$Jaccard_similarity,
        MARGIN = 1,
        FUN = function(x) {
            any(x > 0.99)
        })
    diag(similarities$Jaccard_similarity) <- 1

    similarities$Jaccard_similarity <-
        similarities$Jaccard_similarity[keep_vars, keep_vars]

    keep_vars1 <- apply(X = similarities$is_subset,
        MARGIN = 1,
        FUN = function(x) {
            sum(x) > 1
        })
    keep_vars2 <- apply(X = similarities$is_subset,
        MARGIN = 2,
        FUN = function(x) {
            sum(x) > 1
        })
    keep_vars <- keep_vars1 | keep_vars2

    similarities$is_subset <-
        similarities$is_subset[keep_vars, keep_vars]

    keep_vars1 <- apply(X = similarities$is_almost_subset,
        MARGIN = 1,
        FUN = function(x) {
            sum(x) > 1
        })
    keep_vars2 <- apply(X = similarities$is_almost_subset,
        MARGIN = 2,
        FUN = function(x) {
            sum(x) > 1
        })
    keep_vars <- keep_vars1 | keep_vars2

    similarities$is_almost_subset <-
        similarities$is_almost_subset[keep_vars, keep_vars]

    similarities
}


#' Finds if two vectors only differ between one mutation
#'
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#'
#' @return TRUE, if they only differ by one mutation
differ_by_one_or_less <- function(v1, v2) {
    variants_difference <- v1 == v2
    if (sum(variants_difference) %in% c(length(v1) - 1, length(v1))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Finds the Jaccard similarity between two vectors
#'
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#'
#' @return The Jaccard simularity
jaccard_simularity <- function(v1, v2) {
    variants_difference <- v1 == v2
    return(sum(variants_difference) / length(v1))
}

#' Finds if one variant is a subset of another
#'
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#'
#' @return Result, TRUE if v2 is a subset of v1
is_subset <- function(v1, v2) {
    result <- TRUE
    i <- 1
    while (i <= length(v1) && result == TRUE) {
        if (v1[i] == 0 && v2[i] == 1) {
            return(FALSE)
        }
        i <- i + 1
    }
    return(result)
}

#' Finds if one variant is almost a subset of another
#'
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#'
#' @return Result, TRUE if v2 is almost a subset of v1
is_almost_subset <- function(v1, v2) {
    result <- FALSE
    not_subset_count <- 0
    i <- 1
    while (i <= length(v1)) {
        if (v1[i] == 0) {
            if (v2[i] == 1) {
                not_subset_count <- not_subset_count + 1
            }
        }
        i <- i + 1
    }
    #some work needs to be done on this part
    if ((not_subset_count) / i < 0.005) {
        result <- TRUE
    }
    return(result)
}
