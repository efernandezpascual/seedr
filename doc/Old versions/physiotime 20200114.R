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
#' @export

physiotime <- function(d, t, g, pg, x, reps = NULL, groups = NULL, method)
{
  dtx <- data.table(d)
  if(is.null(reps) & is.null(groups)) {
    dtx[order(get(t)), cumulative := cumsum(get(g)), by = c(x)]
  } else {
    if(is.null(groups)) {
      dtx[order(get(t)), cumulative := cumsum(get(g)), by = c(x, reps)]
    } else {
      if(is.null(reps)) {
        dtx[order(get(t)), cumulative := cumsum(get(g)), by = c(x, groups)]
      } else dtx[order(get(t)), cumulative := cumsum(get(g)), by = c(x, reps, groups)]
    }
  }
  dtx[, germination := cumulative / get(pg)]

  if(is.null(groups)) {
    listx <- list("single group" = dtx)
  } else listx <- split(dtx, by = groups, drop = TRUE)

  listz <- lapply(listx, function(dtx) get(method)(dtx, t, g, pg, x, reps, groups))
  class(listz) <- method
  listz
}

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

bradford <- function(d, t, g, pg, x, reps = NULL, groups = NULL) tryCatch(
  {
    d <- d[! germination %in% c(0, 1)] # Buscar manera sin reasignar
    d[, probit := qnorm(germination, 0, 1)]
    theta <- bradtheta(d[, probit], d[, get(x)], 1/d[, get(t)])
    d[, psibg := get(x) - theta / get(t)]

    hlm <- lm(probit ~ psibg, data = d)
    b <- summary(hlm)$coefficients[1, 1]
    m <- summary(hlm)$coefficients[2, 1]

    z <- list()
    if(is.null(groups)) z$groups <- data.table(group = "single group")
    else {
      z$groups <- unique(d[, ..groups])
      z$groups <- z$groups[, lapply(.SD, as.character)]
      z$groups[, group := do.call(paste, Map(function(x, y) paste(x, y, sep = ":"),
                                             names(.SD), .SD))]
    }

    z$data <- d
    colnames(z$data)[colnames(z$data) == x] <- "psi"
    z$models <- hlm
    if(is.null(reps) & is.null(groups)) {
      z$parameters$n.viable <- sum(d[, .(n = max(get(pg))),
                                       by = .(get(x))][, n])
    } else {
      if(is.null(groups)) {
        z$parameters$n.viable <- sum(d[, .(n = max(get(pg))),
                                         by = .(get(x), get(reps))][, n])
      } else {
        if(is.null(reps)) {
          z$parameters$n.viable <- sum(d[, .(n = max(get(pg))),
                                           by = .(get(x), get(groups))][, n])
        } else z$parameters$n.viable <- sum(d[, .(n = max(get(pg))),
                                                by = .(get(x), get(reps), get(groups))][, n])
      }
    }
    z$parameters$n.germinated <- sum(d[, get(g)])
    z$parameters$levels <- d[order(get(x)), .(Levels = unique(get(x)))][, Levels]
    z$parameters$theta <- theta
    z$parameters$psib50 <- -b/m
    z$parameters$sigma <- 1/m
    z$parameters$R2 <- summary(hlm)$r.squared
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

#' @export
print.bradford <- function(d)
{
  for(y in d) {
    cat("Bradford's hydrotime model", "\n",
        "for the group", y$groups$group, "\n",
        "Number of viable seeds in experiment:",
        y$parameters$n.viable, "\n",
        "Number of germinated seeds in experiment:",
        y$parameters$n.germinated, "\n",
        "Water potential levels in experiment:",
        round(y$parameters$levels, 1), "\n",
        "Theta - Hydrotime constant:",
        round(y$parameters$theta, 2), "\n",
        "Psib50 - Base water potential (median):",
        round((y$parameters$psib50), 2), "\n",
        "Sigma of the base water potential:",
        round(y$parameters$sigma, 2), "\n",
        "R2:", round(y$parameters$R2, 2), "\n", "\n")
  }
}

#' @export
summary.bradford <- function(d)
{
  z <- sapply(d, function(y) {
    data.frame(if(ncol(y$groups) > 1) {
      group = y$groups[,-ncol(y$groups), with = FALSE]
    } else group = y$groups,
    n.viable = y$parameters$n.viable,
    n.germinated = y$parameters$n.germinated,
    n.treatments = length(y$parameters$levels),
    theta = y$parameters$theta,
    psib50 = y$parameters$psib50,
    sigma = y$parameters$sigma,
    R2 = y$parameters$R2)
  })
  dfz <- data.frame(t(z))
  rownames(dfz) <- NULL
  dfz
}

#' @export
plot.bradford <- function(d, treatment) # Plots Bradford's model
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in d) {
    plot(y$data$psibg, y$data$probit, col = as.factor(y$data$psi), pch = 16,
         xlab = expression(paste(psi[b], "(g)")), ylab = "Probit germination")
    abline(lm(y$data$probit ~ y$data$psibg))
    legend("topleft", title = expression(paste(psi, " (MPa)")),
           legend = levels(as.factor(round(y$data$psi, 1))), pch = 16,
           col = unique(as.factor(round(y$data$psi, 1))))
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}
