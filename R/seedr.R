#' seedr: hydro and thermal time seed germination models in R
#'
#' The \code{seedr} package provides functions to fit hydro and thermal time
#' germination models. These models characterize seed lots by two sets of
#' parameters: (i) the physiological thresholds (water, temperature) between
#' which the seed lot can germinate, and (ii) the physiological-time units that
#' the seed lot needs to accumulate before it can germinate. \code{seedr} allows
#' to fit the hydro time model of Bradford (Gummerson 1986, Bradford 1990,
#' Bewley et al. 2013) and the thermal time model of Garcia-Huidobro
#' (Garcia-Huidobro et al. 1982, Gummerson 1986, Bewley et al. 2013).
#' \code{seedr} also allows to quickly fit models to multi-seedlot or
#' multi-species datasets.
#' @references Bewley, J. D., Bradford, K. J., Hilhorst, H. W., & Nonogaki, H.
#'   (2013). Thermal Time Models. In Seeds: Physiology of Development,
#'   Germination and Dormancy, 3rd Edition (pp. 312-317). Springer, New York,
#'   NY.
#'
#'   Bradford, K. J. (1990). A water relations analysis of seed germination
#'   rates. Plant Physiology, 94(2), 840-849.
#'
#'   Garcia-Huidobro, J., Monteith, J. L., & Squire, G. R. (1982). Time,
#'   temperature and germination of pearl millet (Pennisetum typhoides S. & H.)
#'   I. Constant temperature. Journal of Experimental Botany, 33(2), 288-296.
#'
#'   Gummerson, R. J. (1986). The effect of constant temperatures and osmotic
#'   potentials on the germination of sugar beet. Journal of Experimental
#'   Botany, 37(6), 729-741.
#' @import data.table
#' @import graphics
#' @importFrom grDevices colorRampPalette
#' @importFrom stats lm qnorm cov var
#' @docType package
#' @name seedr
NULL

#' Water potential example dataset
#'
#' This is a dataset containing information from a water potential experiment
#' with grass seeds. It is used to give examples of the functions dealing with
#' hydrotime germination models. It also gives and idea of the format in which
#' germination data should be provided to seedr.
#'
#' @format A data frame with 1605 rows and 7 variables \describe{
#'   \item{species}{Name or code for the species to which the data refers}
#'   \item{temperature}{Temperature treatment (in ºC) of the experiment}
#'   \item{psi}{Water potential treatment (in MPa) of the experiment}
#'   \item{dish}{Code for the Petri dish, container, replicate, etc.}
#'   \item{times}{Time of germination scoring since the start of the experiment,
#'   in days} \item{germinated}{Number of germinated seeds recorded at that
#'   time} \item{germinable}{Total number of viables seeds in the replicate} }
#' @source Own data from a laboratory experiment.
"grasses"

#' Temperature example dataset
#'
#' This is a dataset containing information from a germination temperature experiment
#' with centaury seeds. It is used to give examples of the functions dealing with
#' thermal time germination models. It also gives and idea of the format in which
#' germination data should be provided to seedr.
#'
#' @format A data frame with 896 rows and 7 variables \describe{
#'   \item{species}{Name or code for the species to which the data refers}
#'   \item{population}{Name or code for the seedlot}
#'   \item{temperature}{Temperature treatment (in ºC) of the experiment}
#'   \item{dish}{Code for the Petri dish, container, replicate, etc.}
#'   \item{times}{Time of germination scoring since the start of the experiment,
#'   in days} \item{germinated}{Number of germinated seeds recorded at that
#'   time} \item{germinable}{Total number of viables seeds in the replicate} }
#' @source Own data from a laboratory experiment.
"centaury"
