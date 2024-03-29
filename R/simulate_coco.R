#' Simulate counts and coverage from a variant matrix
#' 
#' @param varmat The variant matrix, with properly named rows (VOCs) and columns (mutations)
#' @param rel_counts (Optional) The relative counts of each VOC (will be censored by coverage). Must be same length as nrow(varmat).
#' @param censoring (Optional) The proportion of observations for each mutation. Must be same length as ncol(varmat).
#' @param verbose Print information about the simulation
#' @param absurd If true, the counts are completely unrelated to the variant proportions. Useful for Monte Carlo estimation.
#' 
#' @return A data frame with counts and coverage of the sampled mutations, ready to be used in \code{copt_binom()}.
#' @export
#' 
#' @details By default, uses a negative binomial distribution for the counts and uniform censoring for the coverage.
simulate_coco <- function(varmat, rel_counts = NULL, censoring = NULL, absurd = FALSE, verbose = TRUE) {
    if(is.null(rel_counts)) {
        rel_counts <- stats::rnbinom(nrow(varmat), mu = 25, size = 100)
        max_count <- max(rel_counts)
        rel_counts[sample(1:length(rel_counts), 3, TRUE)] <- max_count * c(3,7,13)
    } else {
        if(length(rel_counts) != nrow(varmat)) {
            stop("Counts must correspond to rows in varmat.")
        }
    }
    if(verbose) {
        print("Expected results:")
        print(rbind(variant = rownames(varmat), 
            proportion = round(rel_counts/sum(rel_counts), 3)))
    }
    if(is.null(censoring)) {
        censoring <- stats::runif(ncol(varmat))
    } else {
        if(length(censoring) != ncol(varmat)) {
            stop("Censoring must correspond to columns in varmat.")
        }
    }

    total_seqs <- sum(rel_counts)
    rho <- rel_counts / total_seqs
    p <- rho %*% varmat
    coverage <- round(censoring * sum(rel_counts))
    counts <- stats::rbinom(n = length(p), size = coverage, prob = p)

    return(data.frame(count = counts, coverage = coverage, mutation = colnames(varmat)))
}
