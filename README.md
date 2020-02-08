
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimRVPedigree

Simulate Pedigrees Ascertained for Rare Disease

# Under Construction

Please note: we are currently expanding `SimRVPedigree` to simulate
diseases with multiple subtypes. For the latest stable version of
`SimRVPedigree` please visit
<https://CRAN.R-project.org/package=SimRVPedigree>.

## Installation

You can install SimRVPedigree from github with:

``` r
# install.packages("devtools")
devtools::install_github("cnieuwoudt/SimRVPedigree")
```

## Overview

`SimRVPedigree` provides routines to simulate and manipulate pedigrees
ascertained to contain multiple family members affected by a rare
disease.

To simulate pedigrees or to simulate life events for an individual we
provide: \* `sim_RVped()` to simulate pedigrees ascertained for multiple
disease-affected relatives, \* `sim_ped()` to simulate pedigrees, and \*
`sim_life()` to simulate life events until a specified stop-year.

To manipulate pedigrees we provide: \* `plot.ped()` to plot pedigrees,
\* `summary.ped()` to obtain summary information for a pedigree or a
collection of pedigrees, \* `censor_ped()` to censor a pedigree after a
specified year, and \* `reassign_gen()` to reassign generation number
based on the most recent common ancestor of all disease-affected
relatives.

## Example

``` r
library(SimRVPedigree)

#Load example age-specific hazard rates data.
data(AgeSpecific_Hazards)

#Create an object of class hazard.
HRates <- hazard(hazardDF = AgeSpecific_Hazards)

#Simulate a pedigree ascertained for multiple disease-affected relatives
set.seed(5444)
ex_ped <- sim_RVped(hazard_rates = HRates,
                    GRR = 20,
                    RVfounder = TRUE,
                    FamID = 1,
                    founder_byears = c(1900, 1950),
                    ascertain_span = c(1995, 2015),
                    num_affected = 2,
                    stop_year = 2017,
                    recall_probs = c(1, 1, 0.25, 0))

# Plot the ascertained pedigree
plot(ex_ped[[2]], ref_year = 2017, cex = 0.75)
```

![](README-example-1.png)<!-- -->

``` r

# Obtain summary information for the ascertained pedigree.
summary(ex_ped[[2]])
#> $family_info
#>   FamID totalRelatives numAffected aveOnsetAge aveIBD ascertainYear segRV
#> 1     1             15           2          56    0.5          2015  TRUE
#> 
#> $affected_info
#>   FamID ID birthYr onsetYr deathYr proband RVstatus
#> 1     1  4    1936    1996    1998   FALSE        1
#> 2     1  7    1963    2015      NA    TRUE        1
```
