#' Simulate counts and coverage from a variant matrix
#' 
#' @param varmat The variant matrix, with properly named rows (VOCs) and columns (mutations)
#' @param rel_counts (Optional) The relative counts of each VOC (will be censored by coverage). Must be same length as nrow(varmat).
#' @param censoring (Optional) The proportion of observations for each mutation. Must be same length as ncol(varmat).
#' 
#' @return A data frame with counts and coverage of the sampled mutations, ready to be used in \code{copt_binom()}.
#' @export
#' 
#' @details By default, uses a negative binomial distribution for the counts and uniform censoring for the coverage.
simulate_coco <- function(varmat, rel_counts = NULL, censoring = NULL) {
    if(is.null(rel_counts)) {
        rel_counts <- stats::rnbinom(nrow(varmat), mu = 25, size = 100)
    } else {
        if(length(rel_counts) != nrow(varmat)) {
            stop("Counts must correspond to rows in varmat.")
        }
    }
    if(is.null(censoring)) {
        censoring <- stats::runif(ncol(varmat))
    } else {
        if(length(censoring) != ncol(varmat)) {
            stop("Censoring must correspond to columns in varmat.")
        }
    }
    # Total number of sequences
    total_seqs <- sum(rel_counts)
    rel_props <- rel_counts / total_seqs

    # Counts of each mutations (after RNA degredation)
    counts <- round(as.numeric(
        t(varmat) %*% rel_counts) * censoring)
    # Coverage (censored by RNA degredation)
    coverage <- round(rep(total_seqs, 
        length(counts)) * censoring)
    return(data.frame(count = counts, coverage = coverage, mutation = colnames(varmat)))
}
