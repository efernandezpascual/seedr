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
