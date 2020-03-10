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

physiotime <- function(d, t = "times", g = "germinated", pg = "germinable", x = "treatment", groups = NULL, method = "bradford")
{
  d <- physiodata(d, t, g, pg, x, groups)
  if(is.null(groups)) {
    z <- get(method)(d$proportions)
  } else {
    listd <- split(d$proportions, by = groups, drop = TRUE)
    z <- lapply(listd, function(d) {
      zi <- get(method)(d)
      zi$groups <- unique(d[, ..groups])
      zi$groups <- zi$groups[, group := do.call(paste, Map(function(x, y)
        paste(x, y, sep = ": "), names(.SD), .SD))]
      zi
    })
    class(z) <- "physiotime"}
  z
}

# physiotime generic functions

#' @export

print.physiotime <- function(d)
{
  cat("A list of physiological time germination models", "\n",
      "calculated for the following", length(d), "groups:", "\n", "\n")
  for(y in d) {
    cat(y$groups$group, "\n")
    print(y)
  }
}

#' @export

summary.physiotime <- function(d)
{
  l <- lapply(d, function(y) data.frame(
    y$groups[,-ncol(y$groups), with = FALSE], summary(y)))
  rbindlist(l)
}

#' @export

plot.physiotime <- function(d, treatment) # Plots Bradford's model
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in d) {
    plot(y)
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}
