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
    } else if (length(shared) / pre < 0.5 && verbose) {
        warning(paste0("Less than ", length(shared) / pre,
            "% of coco's mutations are being used. Consider a larger variant matrix."))
    }
    if (verbose) {
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
variants_similarity <- function(data) {
    
    data <- as.data.frame(t(data))

    
    subset_of_variants <- data |> dplyr::select_if(~ all(. %in% c(0,1)))
    
    similarities <- list()
    
    similarities$Differ_by_one_or_less <- outer(colnames(subset_of_variants), colnames(subset_of_variants), function(x,y) mapply(FUN = differ_by_one_or_less, v1 = subset_of_variants[,x], v2 = subset_of_variants[,y]))
    colnames(similarities$Differ_by_one_or_less) <- colnames(subset_of_variants)
    rownames(similarities$Differ_by_one_or_less) <- colnames(subset_of_variants)
    diag(similarities$Differ_by_one_or_less) <- rep("place holder", nrow(similarities$Differ_by_one_or_less))
    i <- 1
    while (!is.null(nrow(similarities$Differ_by_one_or_less)) && i <= nrow(similarities$Differ_by_one_or_less)) {
      if (sum(similarities$Differ_by_one_or_less[i, ] == rep(FALSE, ncol(similarities$Differ_by_one_or_less))) == ncol(similarities$Differ_by_one_or_less) - 1) {
        similarities$Differ_by_one_or_less <- similarities$Differ_by_one_or_less[-i, , drop = F]
      }
      else {
        i <- i + 1
      }
    }
    if (!("place holder" %in% similarities$Differ_by_one_or_less)) {
      similarities$Differ_by_one_or_less <- NULL
    }
    else{
      i <- 1
      while (!is.null(ncol(similarities$Differ_by_one_or_less)) && i <= ncol(similarities$Differ_by_one_or_less)) {
        if (sum(similarities$Differ_by_one_or_less[,i] == rep(TRUE, nrow(similarities$Differ_by_one_or_less))) == 0) {
          similarities$Differ_by_one_or_less <- similarities$Differ_by_one_or_less[,-i, drop = F]
        } else {
          i <- i + 1
        }
      }
    }
    
    similarities$Differ_by_one_or_less[similarities$Differ_by_one_or_less == "place holder"] <- TRUE
    
    similarities$Jaccard_similarity <- outer(colnames(subset_of_variants), colnames(subset_of_variants), function(x,y) mapply(FUN = jaccard_simularity, v1 = subset_of_variants[,x], v2 = subset_of_variants[,y]))
    colnames(similarities$Jaccard_similarity) <- colnames(subset_of_variants)
    rownames(similarities$Jaccard_similarity) <- colnames(subset_of_variants)
    
    similarities$is_subset <- outer(colnames(subset_of_variants), colnames(subset_of_variants), function(x,y) mapply(FUN = is_subset, v1 = subset_of_variants[,x], v2 = subset_of_variants[,y]))
    colnames(similarities$is_subset) <- colnames(subset_of_variants)
    rownames(similarities$is_subset) <- colnames(subset_of_variants)
    diag(similarities$is_subset) <- rep("place holder", nrow(similarities$is_subset))
    i <- 1
    while (!is.null(nrow(similarities$is_subset)) && i <= nrow(similarities$is_subset)) {
      if (sum(similarities$is_subset[i, ] == rep(FALSE, ncol(similarities$is_subset))) == ncol(similarities$is_subset) - 1) {
        similarities$is_subset <- similarities$is_subset[-i, , drop = F]
      }
      else {
        i <- i + 1
      }
    }
    if (!("place holder" %in% similarities$is_subset)) {
      similarities$is_subset <- NULL
    }
    else{
      i <- 1
      while (!is.null(ncol(similarities$is_subset)) && i <= ncol(similarities$is_subset)) {
        if (sum(similarities$is_subset[,i] == rep(TRUE, nrow(similarities$is_subset))) == 0) {
          similarities$is_subset <- similarities$is_subset[,-i, drop = F]
        } else {
          i <- i + 1
        }
      }
    }
    similarities$is_subset[similarities$is_subset == "place holder"] <- TRUE
    
    similarities$is_almost_subset <- outer(colnames(subset_of_variants), colnames(subset_of_variants), function(x,y) mapply(FUN = is_almost_subset, v1 = subset_of_variants[,x], v2 = subset_of_variants[,y]))
    colnames(similarities$is_almost_subset) <- colnames(subset_of_variants)
    rownames(similarities$is_almost_subset) <- colnames(subset_of_variants)
    diag(similarities$is_almost_subset) <- rep("place holder", nrow(similarities$is_almost_subset))
    i <- 1
    while (!is.null(nrow(similarities$is_almost_subset)) && i <= nrow(similarities$is_almost_subset)) {
      if (sum(similarities$is_almost_subset[i, ] == rep(FALSE, ncol(similarities$is_almost_subset))) == ncol(similarities$is_almost_subset) - 1) {
        similarities$is_almost_subset <- similarities$is_almost_subset[-i, , drop = F]
      }
      else {
        i <- i + 1
      }
    }
    if (!("place holder" %in% similarities$is_almost_subset)) {
      similarities$is_almost_subset <- NULL
    }
    else{
      i <- 1
      while (!is.null(ncol(similarities$is_almost_subset)) && i <= ncol(similarities$is_almost_subset)) {
        if (sum(similarities$is_almost_subset[,i] == rep(TRUE, nrow(similarities$is_almost_subset))) == 0) {
          similarities$is_almost_subset <- similarities$is_almost_subset[,-i, drop = F]
        } else {
          i <- i + 1
        }
      }
    }
    similarities$is_almost_subset[similarities$is_almost_subset == "place holder"] <- TRUE
    
    return(similarities)
}

#' Finds if two vectors only differ between one mutation
#' 
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#' 
#' @return TRUE, if they only differ by one mutation
differ_by_one_or_less <- function(v1,v2) {
  variants_difference <- v1 == v2
  if (sum(variants_difference) == length(v1) - 1 | sum(variants_difference) == length(v1)) {
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#' Finds the Jaccard similarity between two vectors
#' 
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#' 
#' @return The Jaccard simularity
jaccard_simularity <- function(v1,v2) {
  variants_difference <- v1 == v2
  return(sum(variants_difference)/length(v1))
}

#' Finds if one variant is a subset of another
#' 
#' @param v1 vector for comparison
#' @param v2 vector for comparison
#' 
#' @return Result, TRUE if v2 is a subset of v1
is_subset <- function(v1,v2){
  result <- TRUE
  i <- 1
  while (i <= length(v1) & result == TRUE) {
    if (v1[i] == 0) {
      if (v2[i] == 1) {
        return(FALSE)
      }
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
is_almost_subset <- function(v1,v2){
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
  if ((not_subset_count)/i < 0.005) {
    result <- TRUE
  }
  return(result)
}
