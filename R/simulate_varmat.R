#' Simulate a small lineage matrix for testing
#' 
#' Creates a matrix with the specified VoCs and mutation names, but randomly applies mutations. A future version of this function may be useful for simulating mutations in a more meaningful way, but for now it's only useful for testing purposes.
#' 
#' @param VOC Names of lineages
#' @param mutations Names of mutations
#' 
#' @return A properly structured lineage matrix (with rownames).
#' @export
simulate_lineage_defs <- function(lins = c("BA.1", "BA.2", "B.1.1.529"),
    mutations = c("m3037T", "m22599A", "d221943", "i22205GAGCCAGAA", 
        "m24469A", "d219879", "m241T")) {
    lineage_defs <- matrix(
        data = stats::rbinom(n = length(lins) * length(mutations),
            size = 1, prob = 0.7),
        ncol = length(mutations),
        nrow = length(lins)
    )
    colnames(lineage_defs) <- mutations
    rownames(lineage_defs) <- lins
    return(lineage_defs)
}
