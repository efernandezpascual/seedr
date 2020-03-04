# Main input function =======================================

#' Fitting Physiological Time Seed Germination Models.
#'
#' \code{physiotime} fits physiological time models (thermal time, hydrotime) to
#' seed germination data.
#'
#' This function allows to specify the method to be used (e.g. Bradford's
#' hydrotime model). \code{physiotime} is a constructor function that will
#' return a S3 object whose class depends on the method defined to fit the
#' models.
#'
#' @usage physiotime(d, t, g, pg, x, reps = NULL, groups = NULL,
#'   method)
#' @param d a data frame containing the results of a germination experiment.
#'   The data frame should include columns with scoring times, germination
#'   counts (not cumulative), number of potentially germinable seeds, and the
#'   environmental variable of interest. (e.g. temperature or water potential)
#'   (see \code{\link{peg}} example dataset for appropiate structure)
#' @param t the name of a column (quoted) in \code{d} containing a vector of
#'   numeric scoring times
#' @param g the name of a column (quoted) in \code{d} containing a vector of
#'   integer germination counts (non cumulative)
#' @param pg the name of a column (quoted) in \code{d} containing a vector of
#'   integer numbers of potentially germinable seeds
#' @param x the name of a column (quoted) in \code{d} containing a vector of
#'   numeric values for the environmental variable of interest (e.g.
#'   temperature, water potential)
#' @param reps optional, the name of a column (quoted) in \code{d} containing a
#'   vector of replicate information (e.g. Petri dish)
#' @param groups optional, the names of columns (quoted) in \code{d} containing
#'   grouping variables for the experiment that have to be analysed separately
#'   (e.g. different species or populations, different tempperatures in a water
#'   potential experiment, different treatments to break seed dormancy)
#' @param method the method to be used to fit the models, e.g. Bradford's
#'   hydrotime model
#' @return \code{physiotime} returns a S3 object of class equal to the method
#'   specified by the \code{method} argument. The object is a list. The elements
#'   in the first level are the experimental groups as specified by the
#'   \code{groups} argument (if \code{groups = NULL}, there will be a single
#'   first-level element). The second level of the list contains, for each
#'   group, the following elements:
#'   \describe{
#'   \item{groups}{a data.table with
#'   the grouping variables defined by the \code{groups} argument}
#'   \item{data}{a data.table with the original germination data plus new columns with
#'   variables calculated for model fitting}
#'   \item{models}{the model objects
#'   fitted to calculate the physiological time parameters}
#'   \item{parameters}{a list with a summary of the experiment and the physiological time
#'   parameters}}
#' @examples physiotime(d = peg, t = "times", g = "germinated",
#'                      pg = "germinable", x = "psi", reps = "dish",
#'                      groups = c("species", "temperature"),
#'                      method = "bradford")
#' @import data.table
#' @export

physiotime <- function(d, t, g, pg, x, reps = NULL, groups = NULL, method)
{
  d <- physiotable(d, t, g, pg, x, reps, groups)
  if(is.null(groups)) {
    z <- get(method)(d, t, g, pg, x, reps)
  } else {
    listd <- split(d, by = groups, drop = TRUE) # FASTER DATA.TABLE WAY?
    z <- lapply(listd, function(d) {
      zi <- get(method)(d, t, g, pg, x, reps)
      zi$groups <- unique(d[, ..groups])
      zi$groups <- zi$groups[, group := do.call(paste, Map(function(x, y)
        paste(x, y, sep = ": "), names(.SD), .SD))]
      zi
      })
    class(z) <- "physiolist"}
  z
}

# Bradford's hydrotime model =======================================

#' Fitting Bradford's hydrotime model.
#'
#' \code{bradford} is used by \code{\link{physiotime}} when \code{method =
#' "bradford"} to fit a hydrotime seed germination model using the method of
#' Bradford.
#'
#' @usage bradford(d, t, g, pg, x, reps = NULL, groups = NULL)
#' @param d a data.table produced by \code{physiotime} containing the
#'   germination data
#' @inheritParams physiotime
#' @return \code{bradford} returns a list with the results which is passed back
#'   to \code{\link{physiotime}}. The list includes the following elements:
#'   \describe{ \item{groups}{a data.table with the grouping variables defined
#'   by the \code{groups} argument} of \code{physiotime} \item{data}{a
#'   data.table with the original germination data plus new columns with
#'   variables calculated for model fitting} \item{models}{the \code{lm} objects
#'   fitted to calculate the physiological time parameters} \item{parameters}{a
#'   list with a summary of the experiment and the physiological time
#'   parameters}}
#' @export

bradford <- function(d, t, g, pg, x, reps = NULL) tryCatch(
  {
    d <- d[! mean %in% c(0, 1)] # Buscar manera sin reasignar
    d[, probit := qnorm(mean, 0, 1)]
    theta <- bradtheta(d[, probit], d[, get(x)], 1/d[, get(t)])
    d[, psibg := get(x) - theta / get(t)]

    hlm <- lm(probit ~ psibg, data = d)
    b <- summary(hlm)$coefficients[1, 1]
    m <- summary(hlm)$coefficients[2, 1]

    z <- list()

    z$data <- d
    z$proportions <- proportions(d, x, reps)

    z$parameters$n.viable <- sum(z$proportions$germinable)
    z$parameters$n.germinated <- sum(z$proportions$germinated)
    z$parameters$levels <- d[, .(Levels = unique(get(x)))][, Levels]
    z$parameters$theta <- theta
    z$parameters$psib50 <- -b/m
    z$parameters$sigma <- 1/m
    z$parameters$R2 <- summary(hlm)$r.squared

    class(z) <- "bradford"
    z
  }, error = function(e) e)

#' Calculates Bradford's hydrotime constant (Theta).
#'
#' \code{bradtheta} calculates the population hydrotime constant (Theta)
#'   using Bradford's method by solving EXPLICAR LOS CALCULOS de GIL. This
#'   function is mainly used by \code{\link{bradford}} to perform its
#'   calculations.
#'
#' @usage bradtheta(probit, x, invt)
#' @param probit a vector of probit transformations of the cumulative germination proportion
#' @param x a vector of water potentials of the experimental treatment
#' @param invt a vector of the inverse of the germination scoring time since the start of
#'   the experiment
#' @return \code{bradtheta} returns the value of Theta and passes it back to \code{\link{bradford}}
#' @export

bradtheta <- function(probit, x, invt) # Calculates Bradford's Theta
{
  finite <- is.finite(probit) & is.finite(x) & is.finite(invt)
  probit <- probit[finite]
  x <- x[finite]
  invt <- invt[finite]

  Sxy <- cov(probit, x)
  Sxz <- cov(probit, invt)
  Sxx <- var(probit)
  Syz <- cov(x, invt)
  Syy <- var(x)
  Szz <- var(invt)

  a <- -2 * (Sxz ^ 2) * Syz + 2 * Sxy * Sxz * Szz
  b <- 2 * (Sxz ^ 2) * Syy + 4 * Sxy * Syz * Sxz - 2 * (Sxy ^ 2) *
    Szz - 4 * Syz * Sxy * Sxz
  c <- 2 * Syz * Sxy ^ 2 - 2 * Sxy * Sxz * Syy

  max(c(( -b + sqrt(b ^ 2 - 4 * a * c))/(2 * a),
        (-b -sqrt(b ^ 2 - 4 * a * c))/(2 * a)))
}
