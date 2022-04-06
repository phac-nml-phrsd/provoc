# ProVoC

PROportions of Variants of Concern using counts, coverage, and a variant matrix.

- Counts: The number of times a given mutation was observed.
- Coverage: The number of times the position of a given variant was read.
- Variant Matrix: A matrix where the row names are variants, the column names are mutations, and the entries are 1 if the row variant has the column mutation and 0 otherwise.

# Usage

There are two steps to using this software: create the variant matrix and run the model(s).

The variant matrix can be created several ways:

1. User-specified. The columns must be labelled with mutations (in the same format as the ones in your data), and the row names must be the names of the variant.
2. Specify variants of concern, and use Nextstrain Genbank data to determine the mutations and "nuisance" lineages.
    - A nuisance lineage is a lineage that shares a pre-specified percent of mutations with your VoCs.

Data must be specified as follows.
There must be a column labelled `count` and a column labelled `coverage`.
There can be an optional column labelled `sample`, which will be interpreted as a grouping variable and the analysis will be run once for each unique value of `sample`.
Any other columns will be ignored, so you are free to leave them in your data for future purposes (e.g. location and date will still be tied to unique sample ids).

The package currently can only be installed by time travellers who travel forward to a time when I have completed the package.
For all non Time Lords, you can source the `provoc.R` script and use the individual functions in their current state.

Update: Welcome, Time Lords! 
The package now has basic functionality.
Thank you for using your powers/technology to travel half a day into the future.
Next time, I recommend patience.

Here is an overview of the basic functionality of the package.
It may be helpful to look at the structure of the objects created to know what the primary functions expect.

```R
library(provoc)

# Siumlate with defaul parameters
# (Omicron with a small handful of mutations chosen randomly)
varmat <- simulate_varmat()

# Simulate COunts and COverage (with censoring/ RNA degredation)
# Expect 1/6 BA.1, 2/6 BA.2, and 3/6 B.1.1.529
coco <- simulate_coco(varmat, rel_counts = c(100, 200, 300))

copt <- copt_binom(coco, varmat)
copt$par
```

```
## [1] 0.1661582 0.3343673 0.4994745
```

If you have JAGS and `runjags` installed, you can get an estimate of the posterior for each proportion.
The package includes a helper function to put it in a nice format for `ggplot2` faceting.

```R
library(ggplot2)
library(coda) # for gelman.diag()

coda <- coda_binom(coco, varmat)

# coda_binom returns an mcmc.list object.
gelman.diag(coda) 

# provoc includes a nice helper function, 
# which adds `iter` and `chain` columns and (optionally) pivots longer.
codf <- melt_mcmc(coda, pivot = TRUE)

ggplot(codf) + 
    aes(x = iter, y = value, colour = factor(chain)) +
    geom_line() +
    facet_wrap(~ name, scales = "free_y")
```

