# ProVoC

# ProVoC

[![Lifecycle:
development](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental-1)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PROportions of Variants of Concern using counts, coverage, and a variant
matrix.

Builds and diagnoses a model based on:

- Counts: The number of times a given mutation was observed.
- Coverage: The number of times the position of a given variant was
  read.
- Mutation names: whatever format you want, as long as they match the
  names in the Variant Matrix.
- Variant Matrix: A matrix where the row names are variants, the column
  names are mutations, and the entries are 1 if the row variant has the
  column mutation and 0 otherwise.
  - For example: Usher Barcodes
  - All current methods accept fractional entries.

# Usage

There are two steps to using this software: create the mutation
definitions and run the model(s).

If not specified, the mutation definitions uses hardcoded definitions
from the
[cov-lineages/constellations](https://github.com/cov-lineages/constellations)
repo, which contains the representative mutations that were identified
by the PANGO team. The mutation definitions must have names that match
what exists in the data.

``` r
library(provoc)
data(Baaijens)
b1 <- Baaijens [Baaijens$sra == Baaijens$sra[1], ]
b1$mutation <- parse_mutations(b1$label)

head(b1[, c("count", "coverage", "mutation", "label")])
```

      count coverage        mutation   label
    1 14458    14818      aa:S:D614G ~23403G
    2 10431    32699         C12025T ~12025T
    3   759     9577         G29266A ~29266A
    4 23329    23690  aa:orf1a:T265I  ~1059T
    5  6935    32631       aa:M:R44S ~26654T
    6 13866    27715 aa:orf1a:L3352F ~10319T

The main fitting function is designed to mimic `glm()`, with formula
notation that emphasizes the connection to binomial GLM models.

``` r
res <- provoc(cbind(count, coverage) ~ B.1.1.7 + B.1.429 + B.1.617.2,
    data = b1,
    verbose = FALSE)
res
```

    Call:  ~ cbind(count, coverage) B.1.1.7 + B.1.429 + B.1.617.2

    All models converged.

    Top 3 variants:
         rho ci_low ci_high   variant
    2  0.453     NA      NA   B.1.429
    1  0.008     NA      NA   B.1.1.7
    3 <0.001     NA      NA B.1.617.2

We have created a class for `provoc` objects with convenient methods.
For example, plotting the results is achieved as follows:

``` r
plot(res)
```

![](README_files/figure-commonmark/one-sample-res-base-1.png)

We use the convention of
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html) as
a function that creates a `ggplot2` plot based on a particular type of
data. This allows for all of `ggplot2`’s fanciness on top of a pre-made
plot.

``` r
library(ggplot2)
autoplot(res) +
    theme_bw() +
    labs(title = "Results for one sample")
```

![](README_files/figure-commonmark/one-sample-res-ggplot-1.png)

# Multiple Samples

``` r
# First two sampls from Baaijens
b2 <- Baaijens [Baaijens$sra %in% unique(Baaijens$sra)[1:50], ]
b2$mutations <- parse_mutations(b2$label)
head(b1[, c("count", "coverage", "mutation", "label", "sra")])
```

      count coverage        mutation   label         sra
    1 14458    14818      aa:S:D614G ~23403G SRR15505102
    2 10431    32699         C12025T ~12025T SRR15505102
    3   759     9577         G29266A ~29266A SRR15505102
    4 23329    23690  aa:orf1a:T265I  ~1059T SRR15505102
    5  6935    32631       aa:M:R44S ~26654T SRR15505102
    6 13866    27715 aa:orf1a:L3352F ~10319T SRR15505102

Note the “by” argument below.

``` r
res <- provoc(cbind(count, coverage) ~ B.1.1.7 + B.1.429 + B.1.617.2,
    data = b2, by = "sra",
    verbose = FALSE)
res
```

    Call:  ~ cbind(count, coverage) B.1.1.7 + B.1.429 + B.1.617.2

    All models converged.

    Top 6 variants:
          rho ci_low ci_high variant       group
    80  0.979     NA      NA B.1.429 SRR15505128
    37  0.898     NA      NA B.1.1.7 SRR15505114
    83  0.852     NA      NA B.1.429 SRR15505129
    143 0.821     NA      NA B.1.429 SRR15505149
    46  0.783     NA      NA B.1.1.7 SRR15505117
    59  0.699     NA      NA B.1.429 SRR15505121

The plotting functions above work as expected.

``` r
plot(res)
```

![](README_files/figure-commonmark/multi-sample-res-base-1.png)

``` r
autoplot(res)
```

![](README_files/figure-commonmark/multi-sample-res-ggplot-1.png)

In addition, we’ve added special functionality for data with a “date”
column:

``` r
theme_set(theme_bw())
res$date <- lubridate::ymd(b2$date[match(res$group, b2$sra)])
autoplot(res, date_col = "date") + scale_x_date()
```

    Scale for x is already present.
    Adding another scale for x, which will replace the existing scale.

![](README_files/figure-commonmark/multi-date-res-ggplot-1.png)
