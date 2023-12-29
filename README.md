# ProVoC

[![Lifecycle:
development](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental-1)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROportions of Variants of Concern using counts, coverage, and a variant matrix.

Builds and diagnoses a model based on:

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
# (Omicron sublineages with a handful of mutations chosen haphazardly)
varmat <- simulate_varmat()
varmat
```

```
##           m3037T m22599A d221943 i22205GAGCCAGAA m24469A d219879 m241T
## BA.1           1       1       0               0       1       1     0
## BA.2           0       1       0               1       1       0     1
## B.1.1.529      0       1       1               1       1       1     1
```


```R
# Simulate COunts and COverage (with censoring/ RNA degredation)
# Expect 1/6 BA.1, 2/6 BA.2, and 3/6 B.1.1.529
coco <- simulate_coco(varmat, rel_counts = c(100, 200, 300))
coco
```

```
##   count coverage        mutation
## 1    79      401          m3037T
## 2    60       60         m22599A
## 3   246      520         d221943
## 4   407      478 i22205GAGCCAGAA
## 5   473      473         m24469A
## 6   114      176         d219879
## 7   434      517           m241T
```

The `fuse` function is basically a left join on `mutation` that sets up column names nicely. In future iterations of this package this behaviour will be changed. 

The `provoc` function expects column names `count`, `coverage`, and `var_*` (which `fuse()` provides). It then fits the model.

```R
fused <- fuse(coco, varmat)
copt <- provoc(fused = fused, method = "optim")
copt
```

```
##         rho ci_low ci_high   variant sample
## 1 0.1673370     NA      NA      BA.1      1
## 2 0.3573876     NA      NA      BA.2      1
## 3 0.4752754     NA      NA B.1.1.529      1
```


I have included static versions of the "constellations" mutation lists; these are mutation lists that were defined by the team the created the "Omicron", "Delta", etc. naming scheme. They are no longer updated with modern variants, so I've just included them in `constellation_lists`. The variant matrix can be defined as:

```R
varmat <- astronomize()
dim(varmat)
```

```
## [1]  37 325
```

The structure is the exact same as the `varmat` above, but includes every variant up until the first few Omicron recombinants.

The code below will generate data from all variants at random, but add a bunch of Omicron. 

```R
varmat <- varmat[!grepl("+", rownames(varmat), fixed = TRUE), ]
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
## ... Other variants not shown
```

### Multiple Samples

The `provoc()` function looks for a column labelled `sample`, and will apply the analysis separately across each unique value.

```R
rel_counts <- rpois(nrow(varmat), 10)
is_omicron <- rownames(varmat) %in% c("BA.1", "BA.2", "B.1.1.529")
rel_counts[is_omicron] <- c(100, 200, 300)
coco1 <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
coco1$sample <- "1"

rel_counts[is_omicron] <- c(200, 200, 200)
coco2 <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
coco2$sample <- "2"

rel_counts[is_omicron] <- c(300, 200, 100)
coco3 <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
coco3$sample <- "3"

rel_counts[is_omicron] <- c(400, 100, 0)
coco4 <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
coco4$sample <- "4"

coco <- rbind(coco1, coco2) |> rbind(coco3) |> rbind(coco4)

fused <- fuse(coco, varmat)

copt <- provoc(fused = fused, method = "optim")
copt[copt$variant %in% c("BA.1", "BA.2", "B.1.1.529"), ]
```

```
##       rho ci_low ci_high   variant sample
## 5   0.113     NA      NA B.1.1.529      1
## 16  0.223     NA      NA      BA.1      1
## 17  0.331     NA      NA      BA.2      1
## 35  0.222     NA      NA B.1.1.529      2
## 46  0.217     NA      NA      BA.1      2
## 47  0.223     NA      NA      BA.2      2
## 65  0.338     NA      NA B.1.1.529      3
## 76  0.215     NA      NA      BA.1      3
## 77  0.109     NA      NA      BA.2      3
## 95  0.494     NA      NA B.1.1.529      4
## 106 0.125     NA      NA      BA.1      4
## 107 0.000     NA      NA      BA.2      4
```

As expected, the estimated values are close to the truths (but a little bit lower since all of the variants were included). The format of the output makes for easy plotting with `ggplot2`.
