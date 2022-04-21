#' Ensure mutations are present in both coco and varmat
#' 
#' Finds the intersection between the mutations present in coco and varmat. Will squash lineages together if the resulting mutation list is too similar (see details).
#' 
#' @param coco A data frame with a column labelled \code{mutation}.
#' @param varmat Rownames are variants, column names are Mutations.
#' @param min_perc A variant must have at least \code{min_perc} of the mutations.
#' @param vebose Print information about mutations that were removed by the fusion.
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
fuse <- function(coco, varmat, min_perc = 0.01, verbose = TRUE) {
    if(any(colnames(coco) %in% paste0("var_", rownames(varmat)))) {
        stop("coco should not contain column names that are names of lineages. Is this object already fused?")
    }

    pre <- nrow(coco)
    shared <- intersect(coco$mutation, colnames(varmat))
    # We can't say anything about mutations not in varmat.
    coco <- coco[!is.na(coco$mutation),]
    coco <- coco[coco$mutation %in% shared, ]

    if(length(shared) < 3) {
        stop("Too few shared mutations. Are they in the same format?")
    } else if(length(shared) <= 10 & ncol(varmat) > 10) {
        warning("Fewer than 10 shared mutations. Results may be very difficult to interpret.")
    } else if(pre/length(shared) < 0.5) {
        warning(paste0(length(shared)/pre, 
            "% of coco's mutations are being used. Consider a larger variant matrix."))
    }
    if(verbose) {
        perc_rm <- 1 - round(length(shared) / pre, 3)
        print(paste0(100 * perc_rm, 
            "% of the rows of coco have been removed."))
        coco_only <- coco$mutation[!coco$mutation %in% shared]
        varmat_only <- colnames(varmat)[!colnames(varmat) %in% shared]
        print("coco-only mutations removed:")
        print(coco_only)
        print("varmat-only mutations removed")
        print(varmat_only)
    }


    # Lineages with too few mutations (less than 10% are 1s) 
    # TODO: Make this a parameter?
    vari2 <- varmat[, shared]
    too_many_zeros <- apply(vari2, 1, sum) <= (ncol(vari2) * min_perc)
    vari2 <- vari2[!too_many_zeros, ]

    # Squash identical lineages
    i <- 0
    while(i < nrow(vari2)) {
        i <- i + 1
        this_row <- vari2[i, ]
        dupes <- apply(vari2, 1, function(x) mean(this_row == x))
        dupes[i] <- 0
        if(any(dupes > 0.99)) {
            rownames(vari2)[i] <- paste(rownames(vari2)[dupes > 0.99], collapse = "|")
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
    if(!is.null(sample)) fused <- fused[fused$sample == sample, ]
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

