# ProVoC

[![Lifecycle:
development](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental-1)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROportions of Variants of Concern using counts, coverage, and a variant matrix.

- Counts: The number of times a given mutation was observed.
- Coverage: The number of times the position of a given variant was read.
- Mutation names: whatever format you want, as long as they match the names in the Variant Matrix.
- Variant Matrix: A matrix where the row names are variants, the column names are mutations, and the entries are 1 if the row variant has the column mutation and 0 otherwise.

# Usage

There are two steps to using this software: create the variant matrix and run the model(s).

The variant matrix can be created several ways:

1. User-specified. The columns must be labelled with mutations (in the same format as the ones in your data), and the row names must be the names of the variant.
2. Specify variants of concern, and use Nextstrain Genbank data to determine the mutations and "nuisance" lineages.
    - A nuisance lineage is a lineage that shares a pre-specified percent of mutations with your VoCs.
    - This is still a work-in-progress. It might be best to store the mutation lists in a separate location and construct matrices from them.

Data must be specified as follows.
There must be a column labelled `count`, a column labelled `coverage`, and a column labelled `mutation`.
There can be an optional column labelled `sample`, which will be interpreted as a grouping variable and the analysis will be run once for each unique value of `sample`.
Any other columns will be ignored, so you are free to leave them in your data for future purposes (e.g. location and date will still be tied to unique sample ids).

Here is an overview of the basic functionality of the package. TODO: Vignettes.
It may be helpful to look at the structure of the simulated objects to know what the primary functions expect.

```R
library(provoc)

# Siumlate with defaul parameters
# (Omicron sublineages with a handful of mutations chosen randomly)
varmat <- simulate_varmat()

# Simulate COunts and COverage (with censoring/ RNA degredation)
# Expect 1/6 BA.1, 2/6 BA.2, and 3/6 B.1.1.529
coco <- simulate_coco(varmat, rel_counts = c(100, 200, 300))

fused <- fuse(coco, varmat)
copt <- provoc(fused, method = "optim")
copt$par
```

```
## [1] 0.1661582 0.3343673 0.4994745
```

If you have a cloned version of [the constellations repo](https://github.com/cov-lineages/constellations), you can create a variant matrix based on their variant definitions.
This matrix will need be "fused" with coco to ensure that the intersection of the mutation lists is used.

```R
varmat <- astronomize()

rel_counts <- rpois(nrow(varmat), 10)
is_omicron <- rownames(varmat) %in% c("BA.1", "BA.2", "B.1.1.529")
rel_counts[is_omicron] <- c(100, 200, 300)
coco <- simulate_coco(varmat, rel_counts = rel_counts)

fused <- fuse(coco, varmat)

copt <- provoc(fused = fused, method = "optim")
copt[copt$variant %in% c("BA.1", "BA.2", "B.1.1.529"), ]
```

```
##          rho ci_low ci_high   variant sample
## 5  0.1184445     NA      NA B.1.1.529      1
## 16 0.2301833     NA      NA      BA.1      1
## 17 0.3531027     NA      NA      BA.2      1
## ...
```

If you have JAGS and `runjags` installed, you can get an estimate of the posterior for each proportion.
The package includes a helper function to put it in a nice format for `ggplot2` faceting.

```R
library(ggplot2)
library(coda) # for gelman.diag()

coda <- provoc(fused = fused, method = "runjags", quiet = 0)

# coda_binom returns an mcmc.list object; all coda methods will apply
#gelman.diag(coda) # currently broken

# provoc includes a nice helper function,
# which adds `iter` and `chain` columns and (optionally) pivots longer.
# codf <- melt_mcmc(coda, pivot = TRUE) #currently unnecessary

ggplot(coda) +
    aes(x = variant, y = rho, ymin = ci_low, ymax = ci_high) +
    geom_point() + geom_errorbar() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

