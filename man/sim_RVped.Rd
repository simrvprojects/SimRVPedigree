% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RVPedSim_Functions.R
\name{sim_RVped}
\alias{sim_RVped}
\title{Simulate a pedigree ascertained to contain multiple disease-affected relatives}
\usage{
sim_RVped(hazard_rates, GRR, num_affected, ascertain_span, FamID,
  founder_byears, stop_year = NULL, recall_probs = NULL,
  carrier_prob = 0.002, RVfounder = FALSE, NB_params = c(2, 4/7),
  fert = 1, first_diagnosis = NULL, sub_criteria = NULL)
}
\arguments{
\item{hazard_rates}{An object of class \code{hazard}, created by \code{\link{hazard}}.}

\item{GRR}{Numeric. The genetic relative-risk of disease, i.e. the relative-risk of disease for individuals who carry at least one copy of the causal variant.  Note: When simulating diseases with multiple subtypes \code{GRR} must contain one entry for each simulated subtype.  See details.}

\item{num_affected}{Numeric vector.  The minimum number of disease-affected relatives required for ascertainment.}

\item{ascertain_span}{Numeric vector of length 2.  The year span of the ascertainment period.  This period represents the range of years during which the proband developed disease and the family would have been ascertained for multiple affected relatives.}

\item{FamID}{Numeric. The family ID to assign to the simulated pedigree.}

\item{founder_byears}{Numeric vector of length 2.  The span of years from which to simulate, uniformly, the birth year for the founder who introduced the rare variant to the pedigree.}

\item{stop_year}{Numeric. The last year of study.  If not supplied, defaults to the current year.}

\item{recall_probs}{Numeric. The proband's recall probabilities for relatives, see details.  If not supplied, the default value of four times kinship coefficient between the proband and the relative is used.}

\item{carrier_prob}{Numeric.  The carrier probability for all causal variants with relative-risk of disease \code{GRR}.  By default, \code{carrier_prob}\code{ = 0.002}}

\item{RVfounder}{Logical.  Indicates if all pedigrees segregate the rare, causal variant.  By default, \code{RVfounder = FALSE} See details.}

\item{NB_params}{Numeric vector of length 2. The size and probability parameters of the negative binomial distribution used to model the number of children per household.  By default, \code{NB_params}\code{ = c(2, 4/7)}, due to the investigation of Kojima and Kelleher (1962).}

\item{fert}{Numeric.  A constant used to rescale the fertility rate after disease-onset. By default, \code{fert = 1}.}

\item{first_diagnosis}{Numeric. The first year that reliable diagnoses can be obtained regarding disease-affection status.  By default, \code{first_diagnosis}\code{ = NULL} so that all diagnoses are considered reliable. See details.}

\item{sub_criteria}{List. Additional subtype criteria required for ascertainment.  The first item in \code{sub_criteria} is expected to be a character string indicating a subtype label and the second is a numeric entry indicating the minimum number of relatives affected by the identified subtype for ascertianment.   By default, \code{sub_criteria = NULL} so that no additional criteria is applied.  See details.}
}
\value{
A list containing the following data frames:

\item{\code{full_ped} }{The full pedigree, prior to proband selection and trimming.}

\item{\code{ascertained_ped} }{The ascertained pedigree, with proband selected and trimmed according to proband recall probability.  See details.}
}
\description{
\code{sim_RVped} simulates a pedigree ascertained to contain multiple affected members, selects a proband, and trims the pedigree to contain only those individuals that are recalled by the proband.
}
\details{
When \code{RV_founder = TRUE}, all simulated pedigrees will segregate a genetic susceptibility variant.  In this scenario, we assume that the variant is rare enough that it has been introduced by one founder, and we begin the simulation of the pedigree with this founder.  Alternatively, when \code{RV_founder = FALSE} we simulate the starting founder's causal variant status with probability \code{carrier_prob}.  When \code{RV_founder = FALSE} pedigrees may not segregate the genetic susceptibility variant.  The default selection is \code{RV_founder = FALSE}.  Additionally, we note that \code{sim_RVpedigree} is intended for rare causal variants; users will recieve a warning if \code{carrier_prob > 0.002}.

We note that when \code{GRR = 1}, pedigrees do not segregate the causal variant regardless of the setting selected for \code{RVfounder}.  When the causal variant is introduced to the pedigree we transmit it from parent to offspring according to Mendel's laws.

When simulating diseases with multiple subtypes \code{GRR} is a numeric list indicating the genetic-relative risk for each subtype specified in the \code{\link{hazard}} object supplied to \code{hazard_rates}.  For example, for a disease with two disease subtypes, if we set \code{GRR = c(20, 1)} individuals who inherit the causal variant are 20 times more likely than non-carriers to develop the first subtype and as likely as non-carriers to develop the second subtype.

We begin simulating the pedigree by generating the year of birth, uniformly, between the years specified in \code{founder_byears} for the starting founder.  Next, we simulate this founder's life events using the \code{\link{sim_life}} function, and censor any events that occur after the study \code{stop_year}.  Possible life events include: reproduction, disease onset, and death. We continue simulating life events for any offspring, censoring events which occur after the study stop year, until the simulation process terminates.  We do not simulate life events for marry-ins, i.e. individuals who mate with either the starting founder or offspring of the starting founder.

We do not model disease remission. Rather, we impose the restriction that individuals may only experience disease onset once, and remain affected from that point on.  If disease onset occurs then we apply the hazard rate for death in the affected population.

\code{sim_RVped} will only return ascertained pedigrees with at least \code{num_affected} affected individuals.  That is, if a simulated pedigree does not contain at least \code{num_affected} affected individuals \code{sim_RVped} will discard the pedigree and simulate another until the condition is met.  We note that even for \code{num_affected = 2}, \code{sim_RVped} can be computationally expensive.  To simulate a pedigree with no proband, and without a minimum number of affected members use \code{\link{sim_ped}} instead of \code{sim_RVped}.

When simulating diseases with multiple subtypes, users may wish to apply additional ascertainment criteria using the \code{sub_criteria} argument. When supplied, this argument allows users to impose numeric subtype-specific ascertainmet criteria. For example, if and \code{sub_criteria = list("HL", 1)} then at least 1 of the \code{num_affected} disease-affected relatives must be affected by subtype "HL" for the pedigree to be asceratained.  We note that the first entry of \code{sub_criteria}, i.e. the subtype label, must match the one of subtype labels in the hazards object supplied to \code{hazard_rates}.  See examples.

Upon simulating a pedigree with \code{num_affected} individuals, \code{sim_RVped} chooses a proband from the set of available candidates.  Candidates for proband selection must have the following qualities:
\enumerate{
  \item experienced disease onset between the years specified by \code{ascertain_span},
  \item if less than \code{num_affected} - 1 individuals experienced disease onset prior to the lower bound of \code{ascertain_span}, a proband is chosen from the affected individuals, such that there were at least \code{num_affected} affected individuals when the pedigree was ascertained through the proband.
}

We allow users to specify the first year that reliable diagnoses can be made using the argument \code{first_diagnosis}.  All subjects who experience disease onset prior to this year are not considered when ascertaining the pedigree for a specific number of disease-affected relatives.  By default, \code{first_diagnosis = NULL} so that all affected relatives, recalled by the proband, are considered when ascertaining the pedigree.

After the proband is selected, the pedigree is trimmed based on the proband's recall probability of his or her relatives.  This option is included to model the possibility that a proband either cannot provide a complete family history or that they explicitly request that certain family members not be contacted.  If \code{recall_probs} is missing, the default values of four times the kinship coefficient, as defined by Thompson, between the proband and his or her relatives are assumed.  This has the effect of retaining all first degree relatives with probability 1, retaining all second degree relatives with probability 0.5, retaining all third degree relatives with probability 0.25, etc.  Alternatively, the user may specify a list of length \eqn{l}, such that the first \eqn{l-1} items represent the respective recall probabilities for relatives of degree \eqn{1, 2, ... , l-1} and the \eqn{l^{th}} item represents the recall probability of a relative of degree \eqn{l} or greater. For example, if \code{recall_probs = c(1, 0.75, 0.5)}, then all first degree relatives (i.e. parents, siblings, and offspring) are retained with probability 1, all second degree relatives (i.e. grandparents, grandchildren, aunts, uncles, nieces and nephews) are retained with probability 0.75, and all other relatives are retained with probability 0.5. To simulate fully ascertained pedigrees, simply specify \code{recall_probs = c(1)}.

In the event that a trimmed pedigree fails the \code{num_affected} condition,  \code{sim_RVped} will discard that pedigree and simulate another until the condition is met.  For this reason, the values specified for \code{recall_probs} affect computation time.
}
\section{See Also}{

\code{\link{sim_ped}}, \code{\link{trim_ped}}, \code{\link{sim_life}}
}

\examples{
#Read in age-specific hazards
data(AgeSpecific_Hazards)

#Simulate pedigree ascertained for multiple affected individuals
set.seed(2)
ex_RVped <- sim_RVped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
                      GRR = 20,
                      RVfounder = TRUE,
                      FamID = 1,
                      founder_byears = c(1900, 1905),
                      ascertain_span = c(1995, 2015),
                      num_affected = 2,
                      stop_year = 2017,
                      recall_probs = c(1, 1, 0))

# Observe: ex_RVped is a list containing two ped objects
summary(ex_RVped)

# The first is the original pedigree prior
# to proband selection and trimming
plot(ex_RVped[[1]])

# The second is the ascertained pedigree which
# has been trimmed based on proband recall
plot(ex_RVped[[2]])
summary(ex_RVped[[2]])


# NOTE: by default, RVfounder = FALSE.
# Under this setting pedigrees segregate a causal
# variant with probability equal to carrier_prob.



#---------------------------------------------------#
# Simulate Pedigrees with Multiple Disease Subtypes #
#---------------------------------------------------#
# Simulating pedigrees with multiple subtypes
# Import subtype-specific hazards rates for
# Hodgkin's lymphoma and non-Hodgkin's lymphoma
data(SubtypeHazards)
head(SubtypeHazards)

my_hazards <- hazard(SubtypeHazards,
                     subtype_ID = c("HL", "NHL"))


# Simulate pedigree ascertained for at least two individuals
# affected by either Hodgkin's lymphoma or non-Hodgkin's lymphoma.
# Set GRR = c(20, 1) so that individuals who carry a causal variant
# are 20 times more likely than non-carriers to develop "HL" but have
# same risk as non-carriers to develop "NHL".
set.seed(45)
ex_RVped <- sim_RVped(hazard_rates = my_hazards,
                      GRR = c(20, 1),
                      RVfounder = TRUE,
                      FamID = 1,
                      founder_byears = c(1900, 1905),
                      ascertain_span = c(1995, 2015),
                      num_affected = 2,
                      stop_year = 2017,
                      recall_probs = c(1, 1, 0))

plot(ex_RVped[[2]], cex = 0.6)

# Note that we can modify the ascertainment criteria so that
# at least 1 of the two disease-affected relatives are affected by
# the "HL" subtype by supplying c("HL", 1) to the sub_criteria
# argument.
set.seed(69)
ex_RVped <- sim_RVped(hazard_rates = my_hazards,
                      GRR = c(20, 1),
                      RVfounder = TRUE,
                      FamID = 1,
                      founder_byears = c(1900, 1905),
                      ascertain_span = c(1995, 2015),
                      num_affected = 2,
                      stop_year = 2017,
                      recall_probs = c(1, 1, 0),
                      sub_criteria = list("HL", 1))

plot(ex_RVped[[2]], cex = 0.6)

}
\references{
Nieuwoudt, Christina and Jones, Samantha J and Brooks-Wilson, Angela and Graham, Jinko (2018). \emph{Simulating Pedigrees Ascertained for Multiple Disease-Affected Relatives}. Source Code for Biology and Medicine, 13:2.

Ken-Ichi Kojima, Therese M. Kelleher. (1962), \emph{Survival of Mutant Genes}. The American Naturalist 96, 329-346.

Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169.
}
