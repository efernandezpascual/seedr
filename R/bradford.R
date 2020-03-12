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

bradford <- function(d) tryCatch(
  {
    d <- d[! germination.mean %in% c(0, 1)]
    d[, probit := qnorm(germination.mean, 0, 1)]
    theta <- bradtheta(d[, probit], d[, treatment], 1/d[, time])
    d[, psibg := treatment - theta / time]

    hlm <- lm(probit ~ psibg, data = d)
    b <- summary(hlm)$coefficients[1, 1]
    m <- summary(hlm)$coefficients[2, 1]

    z <- list()

    z$data <- d

    z$parameters$levels <- d[, .(Levels = unique(treatment))][, Levels]
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

# bradford generic functions

#' @export

print.bradford <- function(y)
{
  cat("Bradford's hydrotime model", "\n",
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

#' @export

summary.bradford <- function(y)
{
  data.table(
    method = class(y),
    n.treatments = length(y$parameters$levels),
    theta = y$parameters$theta,
    psib50 = y$parameters$psib50,
    sigma = y$parameters$sigma,
    R2 = y$parameters$R2)
}

#' @export

plot.bradford <- function(d) # Plots Bradford's model
{
  colnumber <- length(unique(d$data$treatment))
  colramp <- colorRampPalette(c("red", "orange", "yellow",
                                "green", "blue", "violet"))
  d$data$color <- colramp(colnumber)[as.numeric(cut(d$data$treatment, breaks = colnumber))]

  plot(d$data$psibg, d$data$probit, col = d$data$color, pch = 16,
       xlab = expression(paste(psi[b], " (g)")), ylab = "Probit germination")
  abline(lm(d$data$probit ~ d$data$psibg))
  legend("topleft", title = expression(paste(psi, " (MPa)")),
         legend = levels(as.factor(round(d$data$treatment, 1))), pch = 16,
         col = colramp(colnumber))
}
