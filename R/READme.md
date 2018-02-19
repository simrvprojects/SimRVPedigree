# SimRVPedigree

## Overview
`SimRVPedigree` provides routines to simulate and manipulate pedigrees ascertained to contain multiple family members affected by a rare disease.

To simulate pedigrees or to simulate life events for an individual we provide:
  * `sim_RVped()` to simulate pedigrees ascertained for multiple disease-affected relatives,
* `sim_ped()` to simulate pedigrees, and
* `sim_life()` to simulate life events until a specified stop-year.

To manipulate pedigrees we provide:
  * `plot.ped()` to plot pedigrees,
* `summary.ped()` to obtain summary information for a pedigree or a collection of pedigrees,
* `censor_ped()` to censor a pedigree after a specified year, and
* `reassign_gen()` to reassign generation number based on the most recent common ancestor of all disease-affected relatives.

## Usage
```
library(SimRVPedigree)

#Load example age-specific hazard rates data.
data(AgeSpecific_Hazards)

#Create an object of class hazard.
HRates <- hazard(hazardDF = AgeSpecific_Hazards)

#Simulate a pedigree ascertained for multiple disease-affected relatives
set.seed(1960)
ex_ped <- sim_RVped(hazard_rates = HRates,
                    GRR = 10,
                    RVfounder = TRUE,
                    FamID = 1,
                    founder_byears = c(1900, 1920),
                    ascertain_span = c(1995, 2015),
                    num_affected = 2,
                    stop_year = 2017,
                    recall_probs = c(1, 1, 0))


# Plot pedigree prior to proband selection and trimming
plot(ex_ped[[1]])

# Plot the ascertained pedigree
plot(ex_ped[[2]])

# Obtain summary information for the ascertained pedigree.
summary(ex_ped[[2]])
#> $family_info
#>   FamID totalRelatives numAffected aveOnsetAge aveIBD ascertainYear segRV
#> 1     1             10           3    61.66667    0.5          1999  TRUE
#>
#> $affected_info
#>   FamID ID birthYr onsetYr deathYr RR proband RVstatus
#> 1     1  1    1905    1968    1973 10   FALSE        1
#> 2     1  4    1944    2013      NA 10   FALSE        1
#> 3     1  5    1946    1999    2002  1    TRUE        0
```
