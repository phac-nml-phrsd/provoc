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
notation that emphasizes the connection to binomial GLM models. The
`bootstrap_samples` argument is optional, and is used to calculate
confidence intervals and correlation of the parameters.

``` r
res <- provoc(
    formula = count / coverage ~ B.1.1.7 + B.1.429 + B.1.617.2 + 
        B.1.427 + AY.4 + P.1,
    data = b1,
    bootstrap_samples = 100)
res
```

    Call: count/coverage ~ B.1.1.7 + B.1.429 + B.1.617.2 + B.1.427 + AY.4 + 
        P.1

    Mutations in lineage definitions:  325 
    Mutations used in analysis/mutations in data:
    74/772

    All models converged.

    Top 5 variants:
         rho   ci_low  ci_high   variant group
    5  0.515 5.09e-01 5.20e-01      AY.4     1
    2   0.33 3.23e-01 3.36e-01   B.1.429     1
    4  0.124 1.19e-01 1.30e-01   B.1.427     1
    1  0.008 7.09e-03 7.98e-03   B.1.1.7     1
    3 <0.001 7.45e-10 8.13e-05 B.1.617.2     1

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

``` r
plot_resids(res)
```

![](README_files/figure-commonmark/one-sample-resid-plot-1.png)

``` r
summary(res)
```


    Call:
    count/coverage ~ B.1.1.7 + B.1.429 + B.1.617.2 + B.1.427 + AY.4 + 
        P.1

    Deviance Residuals:
        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    -10.3051  -1.9998   0.0000  -0.7594   0.0000   8.1775 

    Mutations in lineage definitions: 325 
    Mutations used in analysis/mutations in data:
    74/772

    Coefficients:
               rho       ci_low      ci_high   variant
    1 7.599439e-03 7.092760e-03 7.983944e-03   B.1.1.7
    2 3.301205e-01 3.229204e-01 3.355859e-01   B.1.429
    3 9.512833e-09 7.446244e-10 8.130071e-05 B.1.617.2
    4 1.242872e-01 1.193353e-01 1.298808e-01   B.1.427
    5 5.150303e-01 5.088075e-01 5.202288e-01      AY.4

    Correlation of coefficients:
                  B.1.1.7     B.1.429   B.1.617.2     B.1.427        AY.4
    B.1.1.7    1.00000000  0.04022470  0.15339528 -0.09719734  0.07089729
    B.1.429    0.04022470  1.00000000 -0.08722878 -0.88296571  0.06095588
    B.1.617.2  0.15339528 -0.08722878  1.00000000  0.02233181  0.04138764
    B.1.427   -0.09719734 -0.88296571  0.02233181  1.00000000 -0.07725833
    AY.4       0.07089729  0.06095588  0.04138764 -0.07725833  1.00000000

``` r
plot_variants(res)
```

![](README_files/figure-commonmark/one-sample-variant-plot-1.png)

# Multiple Samples

``` r
# First two samples from Baaijens
b2 <- Baaijens [Baaijens$sra %in% unique(Baaijens$sra)[1:30], ]
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

Note the “`by`” argument below. This tells `provoc()` to fit the model
separately to each sample defined by the `by` column. Since there are 30
samples, the model takes a while to run, hence why there are only 50
bootstrap samples here.

``` r
system.time(
    res <- provoc(
        formula = count / coverage ~ B.1.1.7 + B.1.429 + B.1.617.2 +
            B.1.427 + AY.4 + P.1,
        data = b2, 
        by = "sra",
        bootstrap_samples = 0)
)
```

    [1] "Attempt 10 of 20."
    [1] "Attempt 20 of 20."

       user  system elapsed 
      1.456   0.015   1.477 

``` r
res
```

    Call: count/coverage ~ B.1.1.7 + B.1.429 + B.1.617.2 + B.1.427 + AY.4 + 
        P.1

    Mutations in lineage definitions:  325 
    Summary of percent of mutations in data used:
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.09585 0.09585 0.09585 0.09585 0.09585 0.09585 

    All models converged.

    Top 6 variants:
          rho ci_low ci_high variant       group avg_spot_len sample_name     bases
    61  0.899     NA      NA B.1.1.7 SRR15505114          302         FX2 385650980
    15  0.846     NA      NA    AY.4 SRR15505104          302         DR2 412520222
    139  0.84     NA      NA B.1.427 SRR15505129          302         FE1 339620744
    114 0.709     NA      NA B.1.427 SRR15505124          302         FK2 433136554
    50  0.695     NA      NA    AY.4 SRR15505111          302         GC1 484282368
    35  0.687     NA      NA    AY.4 SRR15505108          302         GF2 615765014
         bioproject       date
    61  PRJNA741211 2021-04-11
    15  PRJNA741211 2021-01-15
    139 PRJNA741211 2021-03-16
    114 PRJNA741211 2021-03-24
    50  PRJNA741211 2021-04-17
    35  PRJNA741211 2021-04-23

``` r
plot_variants(res)
```

![](README_files/figure-commonmark/multi-sample-variant-plot-1.png)

The plotting functions above work as expected.

``` r
plot(res)
```

![](README_files/figure-commonmark/multi-sample-res-base-1.png)

``` r
plot_resids(res)
```

![](README_files/figure-commonmark/multi-sample-resid-base-1.png)

``` r
autoplot(res)
```

![](README_files/figure-commonmark/multi-sample-res-ggplot-1.png)

In addition, we’ve added special functionality for data with a “date”
column:

``` r
theme_set(theme_bw())
res$date <- lubridate::ymd(res$date)
autoplot(res, date_col = "date")
```

![](README_files/figure-commonmark/multi-date-res-ggplot-1.png)

``` r
plot_resids(res)
```

![](README_files/figure-commonmark/multi-sample-resid-plot-1.png)

## Searching for Different Variants

``` r
library(patchwork)
res_with <- provoc(cbind(count, coverage) ~ B.1.1.7 + B.1.429 + B.1.617.2 + B.1.427,
    data = b2, by = "sra",
    verbose = FALSE, bootstrap_samples = 0)
res_with$date <- lubridate::ymd(res_with$date)

res_without <- provoc(cbind(count, coverage) ~ B.1.429 + B.1.617.2 + B.1.427,
    data = b2, by = "sra",
    verbose = FALSE, bootstrap_samples = 0)
res_without$date <- lubridate::ymd(res_without$date)

g_with <- autoplot(res_with, date_col = "date") +
    scale_fill_manual(values = 2:5)
g_without <- autoplot(res_without, date_col = "date") +
    scale_fill_manual(values = 3:5)

g_with / g_without
```

![](README_files/figure-commonmark/with-without-1.png)

# Roadmap

- Soon
  - [ ] Finish simplified versions of Freyja and Alcov models.
  - [ ] Full suite of model diagnostic visualizations.
  - [ ] Process data from iVar for use in models.
    - [ ] Process a list or directory of iVar data into an appropriate
      object.
  - [ ] Tools to investigate the effect of different lineage
    definitions.
- Medium future
  - [ ] Convert Freyja and Alcov outputs to R objects that can take
    advantage of the diagnostics.
  - [ ] Allow lists of formulae, lists of mutation_defs, and lists of
    data.
- Far future
  - [ ] Make available on CRAN
