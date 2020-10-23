#' seedr: A package to fit hydro and thermal time seed germination models
#'
#' The seedr package provides functions to process data generated in laboratory
#' germination experiments, fitting hydrothermal time germination models. These
#' models allow to characterize seed lots by parameters that indicate the physiological
#' thresholds (of temperature or water) between which the seed lot can germinate,
#' and the physiological-time units that the seeds in the seed lot need to accumulate
#' before they can germinate.
#'
#' @import data.table
#' @docType package
#' @name seedr
NULL

#' Water potential example dataset.
#'
#' This is a dataset containing information from a water potential experiment
#'   with grass seeds. It is used to give examples of the functions dealing with
#'   hydrotime germination models. It also shows the format in which germination
#'   data should be provided to seedr.
#'
#' @format A data frame with 384 rows and 8 variables
#' \describe{
#'   \item{species}{Name or code for the species to which the data refers}
#'   \item{population}{Name or code for the population or seed accession}
#'   \item{temperature}{Temperature treatment (in ºC) of the experiment}
#'   \item{psi}{Water potential treatment (in MPa) of the experiment}
#'   \item{dish}{Code for the Petri dish, container, replicate, etc.}
#'   \item{times}{Time of germination scoring since the start of the experiment,
#'     in days}
#'   \item{germinated}{Number of germinated seeds recorded at that time}
#'   \item{germinable}{Total number of viables seeds in the replicate}
#'   }
#' @source Own data from a laboratory experiment.
"grasses"
