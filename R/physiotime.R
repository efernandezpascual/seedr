#' Fits physiological time seed germination models
#'
#' \code{physiotime} fits physiological time models (thermal time, hydrotime) to
#' seed germination data. It is a wrapper function that transforms data to class
#' "physiodata" and allows to specify the physiological time model to be fitted
#' (i.e. Bradford's hydrotime model or Garcia-Huidobro's thermal time model).
#'
#' @usage physiotime(d, t = "times", g = "germinated", pg = "germinable", x =
#'   "treatment", groups = NULL, method = "bradford", min.ptos = 3, tops =
#'   c("Max R2","Max value"), fractions = (1:9)/10)
#' @param d a data.frame containing the results of a germination experiment. The
#'   data frame should include columns with scoring times, germination counts
#'   (not cumulative), number of potentially germinable seeds, and the
#'   environmental variable of interest. (e.g. temperature or water potential)
#'   (see \code{\link{grasses}} example dataset for appropiate structure).
#' @param t the name of a column in \code{d} containing a vector of numeric
#'   scoring times.
#' @param g the name of a column in \code{d} containing a vector of integer
#'   germination counts (non cumulative).
#' @param pg the name of a column in \code{d} containing a vector of integer
#'   numbers of potentially germinable seeds.
#' @param x the name of a column in \code{d} containing a vector of numeric
#'   values for the environmental variable of interest (e.g. temperature, water
#'   potential).
#' @param groups optional, the names of columns in \code{d} containing grouping
#'   variables for the experiment that have to be analysed separately (e.g.
#'   different species or populations, different temperatures in a water
#'   potential experiment, different treatments to break seed dormancy).
#' @param method the method to be used to fit the models, can be "bradford" to
#'   fit a hydrotime model or "huidobro" to fit a thermal time model.
#' @param min.ptos minimal number of data points (i.e. different temperature
#'   treatments) needed to fit the suboptimal and supraoptimal germination
#'   models if fitting a thermal time model. If the number of points available
#'   in the dataset is less than \code{min.ptos}, then the suboptimal or the
#'   supraoptimal models are not fitted.
#' @param tops method used to divide the dataset in suboptimal and supraoptimal
#'   sections if fitting a thermal time model. "Max value" splits the data by
#'   the temperature that produces the highest seed germination rate. "Max R2"
#'   splits the data by the temperature that maximises the R2 of the suboptimal
#'   and supraoptimal linear regressions.
#' @param fractions percentiles into which the seed population is split if
#'   fitting a thermal time model. The default is the 9 deciles (i.e. t10, t20..
#'   t90) as used by Garcia-Huidobro.
#' @return \code{physiotime} returns a S3 object of class "physiotime".
#'  The object is a list containing, for each group (seedlot, speciec, etc.) the
#'  results of fitting the physiological time models.  The generic functions
#'   \code{summary} and \code{plot} are used to obtain and visualize the model
#'   results.
#' @examples
#' m <- physiotime(centaury, x = "temperature",
#'                 method = "huidobro", groups = c("species", "population"))
#' m
#' summary(m)
#' plot(m)
#' @export
physiotime <- function(d, t = "times", g = "germinated", pg = "germinable",
                       x = "treatment", groups = NULL, method = "bradford",
                       min.ptos = 3, tops = c("Max R2","Max value"),
                       fractions = (1:9)/10)
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
print.physiotime <- function(x, ...)
{
  cat("A list of physiological time germination models", "\n",
      "calculated for the following", length(x), "groups:", "\n", "\n")
  for(y in x) {
    cat(y$groups$group, "\n")
    print(y)
  }
}

#' @export
summary.physiotime <- function(object, ...)
{
  l <- lapply(object, function(y) data.frame(
    y$groups[,-ncol(y$groups), with = FALSE], summary(y)))
  rbindlist(l)
}

#' @export
plot.physiotime <- function(x, ...) # Plots Bradford's model
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in x) {
    plot(y)
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}
