#' Simulate counts and coverage from a lineage definition matrix
#' 
#' @param lineage_defs The lineage definition matrix, with properly named rows (VOCs) and columns (mutations)
#' @param rel_counts (Optional) The relative counts of each VOC (will be censored by coverage). Must be same length as nrow(lineage_defs).
#' @param censoring (Optional) The proportion of observations for each mutation. Must be same length as ncol(lineage_defs).
#' @param verbose Print information about the simulation
#' @param absurd If true, the counts are completely unrelated to the lineage proportions. Useful for Monte Carlo estimation.
#' 
#' @return A data frame with counts and coverage of the sampled mutations, ready to be used in \code{copt_binom()}.
#' @export
#' 
#' @details By default, uses a negative binomial distribution for the counts and uniform censoring for the coverage.
simulate_coco <- function(lineage_defs, rel_counts = NULL, censoring = NULL, absurd = FALSE, verbose = TRUE) {
    if(is.null(rel_counts)) {
        rel_counts <- stats::rnbinom(nrow(lineage_defs), mu = 25, size = 100)
        max_count <- max(rel_counts)
        rel_counts[sample(1:length(rel_counts), 3, TRUE)] <- max_count * c(3,7,13)
    } else {
        if(length(rel_counts) != nrow(lineage_defs)) {
            stop("Counts must correspond to rows in lineage_defs.")
        }
    }
    if(verbose) {
        print("Expected results:")
        print(rbind(lineage = rownames(lineage_defs), 
            proportion = round(rel_counts/sum(rel_counts), 3)))
    }
    if(is.null(censoring)) {
        censoring <- stats::runif(ncol(lineage_defs))
    } else {
        if(length(censoring) != ncol(lineage_defs)) {
            stop("Censoring must correspond to columns in lineage_defs.")
        }
    }

    total_seqs <- sum(rel_counts)
    rho <- rel_counts / total_seqs
    p <- rho %*% as.matrix(lineage_defs)
    coverage <- round(censoring * sum(rel_counts))
    counts <- stats::rbinom(n = length(p), size = coverage, prob = p)

    return(data.frame(count = counts, coverage = coverage, mutation = colnames(lineage_defs)))
}
