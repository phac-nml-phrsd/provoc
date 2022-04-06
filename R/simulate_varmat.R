#' Simulate a small variant matrix for testing
#' 
#' @param VOC Names of variants of concern
#' @param mutations Names of mutations
#' 
#' @return A properly structured variant matrix (with rownames).
#' @export
simulate_varmat <- function(VOC = c("BA.1", "BA.2", "B.1.1.529"),
    mutations = c("m3037T", "m22599A", "d221943", "i22205GAGCCAGAA", 
        "m24469A", "d219879", "m241T")) {
    varmat = matrix(
        data = stats::rbinom(n = length(VOC) * length(mutations), 
            size = 1, prob = 0.7),
        ncol = length(mutations),
        nrow = length(VOC)
    )
    colnames(varmat) <- mutations
    rownames(varmat) <- VOC
    return(varmat)
}
