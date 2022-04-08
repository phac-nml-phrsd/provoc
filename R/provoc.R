#' Proportions of Variants of Concern
#' 
#' Un-fuses coco and varmat and applies the appropriate estimation technique. Applies to each unique value of 
#' 
#' @param fused The fused data frame of coco and varmat
#' @param method Optimization or Bayesian?
#' 
#' @return Results of the estimation. Currently has different output depending on method, but will hopefully be unified in a later version.
#' @export
provoc <- function (fused, method = c("optim", "runjags")) {
    # TODO: Respect the sample column (output a list.)
    # TODO: Allow parameters to be passed to the respective functions.
    variants <- startsWith(names(fused), "var_")
    coco <- fused[, !variants]

    vardf <- fused[, variants]
    varmat <- t(as.matrix(vardf))
    varmat <- matrix(as.numeric(varmat), ncol = ncol(varmat))
    rownames(varmat) <- gsub("var_", "", names(coco)[variants])
    colnames(varmat) <- coco$mutation

    if(method[1] == "optim") {
        res <- provoc_optim(coco, varmat)
        # TODO: Process the output
    } else {
        res <- provoc_jags(coco, varmat)
        # TODO: Process the output
    }
    # TODO: Unify output formats, including convergence yes/no, convergence note and point estimates.
    # TODO: (Long Term) make this it's own class with nice printing defaults
    res
}



