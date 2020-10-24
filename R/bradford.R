#' Fits Bradford's hydrotime model
#'
#' \code{bradford} fits a hydrotime seed germination model using the method of
#' Bradford (Gummerson 1986, Bradford 1990, Bewley et al. 2013). This function
#' can be used only with one-group dataset, i.e. one seed lot of one species. To
#' fit models to grouped datasets (multi-seedlots, multi-species) use the
#' function \code{physiotime} instead.
#'
#' @usage bradford(d)
#' @param d a data.table within a "physiodata" object, containing the cumulative
#'   germination proportion at each scoring time and water potential treatment.
#' @return \code{bradford} returns a S3 object of class "bradford" with the
#'   results of fitting the hydrotime model. The generic functions
#'   \code{summary} and \code{plot} are used to obtain and visualize the model
#'   results.
#' @examples
#' # format dataset with physiodata
#' anisantha <- physiodata(subset(grasses, species == "Anisantha rubens"), x = "psi")
#' # bradford() uses the $proportions element within the physiodata object
#' b <- bradford(anisantha$proportions)
#' b # prints the main hydrotime variables
#' summary(b) # returns the main hydrotime variables as a data.table
#' plot(b) # plots the fitted model
#' @references Bewley, J. D., Bradford, K. J., Hilhorst, H. W., & Nonogaki, H.
#'   (2013). Hydrotime Model of Germination. In Seeds: Physiology of
#'   Development, Germination and Dormancy, 3rd Edition (pp. 303-307). Springer,
#'   New York, NY.
#'
#'   Bradford, K. J. (1990). A water relations analysis of seed germination
#'   rates. Plant Physiology, 94(2), 840-849.
#'
#'   Gummerson, R. J. (1986). The effect of constant temperatures and osmotic
#'   potentials on the germination of sugar beet. Journal of Experimental
#'   Botany, 37(6), 729-741.
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

# bradford generic functions

#' @export
print.bradford <- function(x, ...)
{
  cat("Bradford's hydrotime model", "\n",
      "Water potential levels in experiment:",
      round(x$parameters$levels, 1), "\n",
      "Theta - Hydrotime constant:",
      round(x$parameters$theta, 2), "\n",
      "Psib50 - Base water potential (median):",
      round((x$parameters$psib50), 2), "\n",
      "Sigma of the base water potential:",
      round(x$parameters$sigma, 2), "\n",
      "R2:", round(x$parameters$R2, 2), "\n", "\n")
}

#' @export
summary.bradford <- function(object, ...)
{
  data.table(
    method = class(object),
    n.treatments = length(object$parameters$levels),
    theta = object$parameters$theta,
    psib50 = object$parameters$psib50,
    sigma = object$parameters$sigma,
    R2 = object$parameters$R2)
}

#' @export
plot.bradford <- function(x, ...) # Plots Bradford's model
{
  colnumber <- length(unique(x$data$treatment))
  colramp <- colorRampPalette(c("red", "orange", "yellow",
                                "green", "blue", "violet"))
  x$data$color <- colramp(colnumber)[as.numeric(cut(x$data$treatment, breaks = colnumber))]

  plot(x$data$psibg, x$data$probit, col = x$data$color, pch = 16,
       xlab = expression(paste(psi[b], " (g)")), ylab = "Probit germination")
  abline(lm(x$data$probit ~ x$data$psibg))
  legend("topleft", title = expression(paste(psi, " (MPa)")),
         legend = levels(as.factor(round(x$data$treatment, 1))), pch = 16,
         col = colramp(colnumber))
}

# bradford internal functions

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
