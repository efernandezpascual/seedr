#' Fits Garcia-Huidobro's thermal time model
#'
#' \code{huidobro} fits a thermal time seed germination model using the method
#' of Garcia-Huidobro (Garcia-Huidobro et al. 1982, Gummerson 1986, Bewley et
#' al. 2013). This function can be used only with one-group dataset, i.e. one
#' seed lot of one species. To fit models to grouped datasets (multi-seedlots,
#' multi-species) use the function \code{physiotime} instead.
#'
#' @usage huidobro(d, min.ptos = 3, tops = c("Max R2","Max value"), fractions =
#'   (1:9)/10)
#' @param d a data.table within a "physiodata" object, containing the cumulative
#'   germination proportion at each scoring time and temperature treatment.
#' @param min.ptos minimal number of data points (i.e. different temperature
#'   treatments) needed to fit the suboptimal and supraoptimal germination
#'   models. If the number of points available in the dataset is less than
#'   \code{min.ptos}, then the suboptimal or the supraoptimal models are not
#'   fitted.
#' @param tops method used to divide the dataset in suboptimal and supraoptimal
#'   sections. "Max value" splits the data by the temperature that produces the
#'   highest seed germination rate. "Max R2" splits the data by the temperature
#'   that maximises the R2 of the suboptimal and supraoptimal linear
#'   regressions.
#' @param fractions percentiles into which the seed population is split to fit
#'   the thermal time model. The default is the 9 deciles (i.e. t10, t20.. t90)
#'   as used by Garcia-Huidobro.
#' @return \code{huidobro} returns a S3 object of class "huidobro" with the
#'   results of fitting the thermal time model. The generic functions
#'   \code{summary} and \code{plot} are used to obtain and visualize the model
#'   results.
#' @examples
#' # format dataset with physiodata
#' malva <- physiodata(subset(centaury, population == "La Malva"), x = "temperature")
#' # huidobro() uses the $proportions element within the physiodata object
#' h <- huidobro(malva$proportions)
#' h # prints the main thermal time variables
#' summary(h) # returns the main thermal time variables as a data.table
#' plot(h) # plots the fitted model
#' @references Bewley, J. D., Bradford, K. J., Hilhorst, H. W., & Nonogaki, H.
#'   (2013). Thermal Time Models. In Seeds: Physiology of Development,
#'   Germination and Dormancy, 3rd Edition (pp. 312-317). Springer, New York,
#'   NY. Bradford, K. J. (1990). A water relations analysis of seed germination
#'   rates. Plant Physiology, 94(2), 840-849.
#'
#'   Garcia-Huidobro, J., Monteith, J. L., & Squire, G. R. (1982). Time,
#'   temperature and germination of pearl millet (Pennisetum typhoides S. & H.)
#'   I. Constant temperature. Journal of Experimental Botany, 33(2), 288-296.
#'
#'   Gummerson, R. J. (1986). The effect of constant temperatures and osmotic
#'   potentials on the germination of sugar beet. Journal of Experimental
#'   Botany, 37(6), 729-741.
#' @export
huidobro <- function(d, min.ptos = 3, tops = c("Max R2","Max value"),
                     fractions = (1:9)/10) tryCatch(
                       {
                         rates <- d[, .(fraction = fractions, rate = rates(d = .SD, fractions = fractions)), by = c("treatment")]
                         rates <- rates[! is.na(rates$rate)]

                         if ("Max value" %in% tops)  # Tb<=Tmin, Tc>=Tmax for all quantiles
                         {
                           Ts <- rates[, cardinals(.SD[["treatment"]], rate, which.max(rate),
                                                   min.ptos = min.ptos), by = c("fraction")]
                           Topt <- Ts[, .(Tb = min(c(mean(Tb, na.rm = TRUE), min(Tmin, na.rm = TRUE))),
                                          Tc = max(c(mean(Tc, na.rm = TRUE), max(Tmax), na.rm = TRUE)),
                                          To = mean(To, na.rm = TRUE))]
                         }

                         if ("Max R2" %in% tops)  # Tb<=Tmin, Tc>=Tmax for all quantiles
                         {
                           Ts <- rates[, cardinals(.SD[["treatment"]], rate,
                                                   maximizeR2(.SD[["treatment"]], rate, min.ptos = min.ptos),
                                                   min.ptos = min.ptos), by = c("fraction")]
                           Topt <- Ts[, .(Tb = min(c(mean(Tb, na.rm = TRUE), min(Tmin, na.rm = TRUE))),
                                          Tc =max(c(mean(Tc, na.rm = TRUE), max(Tmax), na.rm = TRUE)),
                                          To = mean(To, na.rm = TRUE))]
                         }

                         d[, probit := qnorm(germination.mean, 0, 1)]

                         subs <- d[treatment <= Topt$To & is.finite(probit)]
                         subs[, thetag := time * (treatment - Topt$Tb)]
                         subfit <- gaussfit(subs)
                         sups <- d[treatment >= Topt$To & is.finite(probit)]
                         sups[, thetag := time * (Topt$Tc - treatment)]
                         supfit <- gaussfit(sups)

                         l <- list()

                         l$data$suboptimal <- subs
                         l$data$supraoptimal <- sups

                         l$parameters$levels$suboptimal <- unique(subs$treatment)
                         l$parameters$levels$supraoptimal <- unique(sups$treatment)
                         l$parameters$cardinals <- Topt
                         l$parameters$suboptimal <- subfit
                         l$parameters$supraoptimal <- supfit

                         class(l) <- "huidobro"
                         l
                       }, error = function(e) e)

# huidobro generic functions

#' @export
print.huidobro <- function(x, ...)
{
  if(is.nan(x$parameters$cardinals$Tc)) {
    cat("Garcia-Huidobro's suboptimal thermal time model", "\n",
        "Suboptimal temperature levels in experiment:",
        round(x$parameters$levels$suboptimal, 1), "\n",
        "Tb - Base temperature:",
        round(x$parameters$cardinals$Tb, 1), "\n",
        "theta50 1 - Suboptimal thermal time (median):",
        round((x$parameters$suboptimal$theta50), 2), "\n",
        "Sigma of the suboptimal thermal time:",
        round((x$parameters$suboptimal$sigma), 2), "\n",
        "R2 of the suboptimal model:",
        round(x$parameters$suboptimal$R2, 2), "\n")} else {
          if(is.nan(x$parameters$cardinals$To)) {
            cat("Garcia-Huidobro's supraoptimal thermal time model", "\n",
                "Supraoptimal temperature levels in experiment:",
                round(x$parameters$levels$supraoptimal, 1), "\n",
                "Tc - Ceiling temperature:",
                round(x$parameters$cardinals$Tc, 1), "\n",
                "theta50 1 - Suboptimal thermal time (median):",
                round(x$parameters$supraoptimal$theta50, 2), "\n",
                "Sigma of the suboptimal thermal time:",
                round(x$parameters$supraoptimal$sigma, 2), "\n",
                "R2 of the suboptimal model:",
                round(x$parameters$supraoptimal$R2, 2), "\n", "\n")} else
                  cat("Garcia-Huidobro's thermal time model", "\n",
                      "Suboptimal temperature levels in experiment:",
                      round(x$parameters$levels$suboptimal, 1), "\n",
                      "Supraoptimal temperature levels in experiment:",
                      round(x$parameters$levels$supraoptimal, 1), "\n",
                      "Tb - Base temperature:",
                      round(x$parameters$cardinals$Tb, 1), "\n",
                      "To - Optimal temperature:",
                      round(x$parameters$cardinals$To, 1), "\n",
                      "Tc - Ceiling temperature:",
                      round(x$parameters$cardinals$Tc, 1), "\n",
                      "theta50 1 - Suboptimal thermal time (median):",
                      round((x$parameters$suboptimal$theta50), 2), "\n",
                      "Sigma of the suboptimal thermal time:",
                      round((x$parameters$suboptimal$sigma), 2), "\n",
                      "R2 of the suboptimal model:",
                      round(x$parameters$suboptimal$R2, 2), "\n",
                      "theta50 2 - Supraoptimal thermal time (median):",
                      round(x$parameters$supraoptimal$theta50, 2), "\n",
                      "Sigma of the supraoptimal thermal time:",
                      round(x$parameters$supraoptimal$sigma, 2), "\n",
                      "R2 of the supraoptimal model:",
                      round(x$parameters$supraoptimal$R2, 2), "\n", "\n")}
}

#' @export
summary.huidobro <- function(object, ...)
{
  data.table(
    method = class(object),
    n.treatments.sub = length(object$parameters$levels$suboptimal),
    n.treatments.sup = length(object$parameters$levels$supraoptimal),
    Tb = object$parameters$cardinals$Tb,
    To = object$parameters$cardinals$To,
    Tc = object$parameters$cardinals$Tc,
    theta50.sub = object$parameters$suboptimal$theta50,
    sigma.sub = object$parameters$suboptimal$sigma,
    R2.sub = object$parameters$suboptimal$R2,
    theta50.sup = object$parameters$supraoptimal$theta50,
    sigma.sup = object$parameters$supraoptimal$sigma,
    R2.sup = object$parameters$supraoptimal$R2)
}

#' @export
plot.huidobro <- function(x, ...) # Plots García-Huidobro's model
{
  if(is.nan(x$parameters$cardinals$Tc)) {
    xpd.status <- par()$xpd
    mar.status <- par()$mar

    par(mar = mar.status + c(0, 0, 0, 4))
    colnumber <- length(unique(x$data$suboptimal$treatment))
    colramp <- colorRampPalette(c("red", "orange", "yellow",
                                  "green", "blue", "violet"))
    x$data$suboptimal$color <- colramp(colnumber)[as.numeric(cut(x$data$suboptimal$treatment, breaks = colnumber))]

    plot(x$data$suboptimal$thetag, x$data$suboptimal$probit, col = x$data$suboptimal$color, pch = 16,
         xlab = expression(paste("Suboptimal ", theta, " (g)")), ylab = "Probit germination")
    abline(lm(x$data$suboptimal$probit ~ x$data$suboptimal$thetag))
    par(xpd = TRUE)
    legend(max(x$data$suboptimal[, thetag]) + max(x$data$suboptimal[, thetag])*.05, 2,
           title = "U+00B0C",
           legend = levels(as.factor(round(x$data$suboptimal$treatment, 1))), pch = 16,
           col = colramp(colnumber))
    par(xpd = xpd.status)
    par(mar = mar.status)} else
    {if(is.nan(x$parameters$cardinals$To)) {
      xpd.status <- par()$xpd
      mar.status <- par()$mar

      par(mar = mar.status + c(0, 0, 0, 4))
      colnumber <- length(unique(x$data$supraoptimal$treatment))
      colramp <- colorRampPalette(c("red", "orange", "yellow",
                                    "green", "blue", "violet"))
      x$data$supraoptimal$color <- colramp(colnumber)[as.numeric(cut(x$data$supraoptimal$treatment, breaks = colnumber))]

      plot(x$data$supraoptimal$thetag, x$data$supraoptimal$probit, col = x$data$supraoptimal$color, pch = 16,
           xlab = expression(paste("Supraoptimal ", theta, " (g)")), ylab = "Probit germination")
      abline(lm(x$data$supraoptimal$probit ~ x$data$supraoptimal$thetag))
      par(xpd = TRUE)
      legend(max(x$data$supraoptimal[, thetag]) + max(x$data$supraoptimal[, thetag])*.05, 2,
             title = "U+00B0C",
             legend = levels(as.factor(round(x$data$supraoptimal$treatment, 1))), pch = 16,
             col = colramp(colnumber))
      par(xpd = xpd.status)
      par(mar = mar.status)} else {
        xpd.status <- par()$xpd
        mar.status <- par()$mar
        ask.status <- par()$ask
        par(ask = TRUE)
        par(mar = mar.status + c(0, 0, 0, 4))

        colnumber <- length(unique(x$data$suboptimal$treatment))
        colramp <- colorRampPalette(c("red", "orange", "yellow",
                                      "green", "blue", "violet"))
        x$data$suboptimal$color <- colramp(colnumber)[as.numeric(cut(x$data$suboptimal$treatment, breaks = colnumber))]

        plot(x$data$suboptimal$thetag, x$data$suboptimal$probit, col = x$data$suboptimal$color, pch = 16,
             xlab = expression(paste("Suboptimal ", theta, " (g)")), ylab = "Probit germination")
        abline(lm(x$data$suboptimal$probit ~ x$data$suboptimal$thetag))
        par(xpd = TRUE)
        legend(max(x$data$suboptimal[, thetag]) + max(x$data$suboptimal[, thetag])*.05, 2,
               title = "U+00B0C",
               legend = levels(as.factor(round(x$data$suboptimal$treatment, 1))), pch = 16,
               col = colramp(colnumber))
        par(xpd = xpd.status)

        colnumber <- length(unique(x$data$supraoptimal$treatment))
        colramp <- colorRampPalette(c("red", "orange", "yellow",
                                      "green", "blue", "violet"))
        x$data$supraoptimal$color <- colramp(colnumber)[as.numeric(cut(x$data$supraoptimal$treatment, breaks = colnumber))]

        plot(x$data$supraoptimal$thetag, x$data$supraoptimal$probit, col = x$data$supraoptimal$color, pch = 16,
             xlab = expression(paste("Supraoptimal ", theta, " (g)")), ylab = "Probit germination")
        abline(lm(x$data$supraoptimal$probit ~ x$data$supraoptimal$thetag))
        par(xpd = TRUE)
        legend(max(x$data$supraoptimal[, thetag]) + max(x$data$supraoptimal[, thetag])*.05, 2,
               title = "U+00B0C",
               legend = levels(as.factor(round(x$data$supraoptimal$treatment, 1))), pch = 16,
               col = colramp(colnumber))

        par(xpd = xpd.status)
        par(mar = mar.status)
        par(ask = ask.status)
      }}
}

# huidobro internal functions

cardinals <- function(x, r, pos, min.ptos = 3)
{
  n <- length(r)
  a1 <- NULL
  b1 <- NULL

  if (pos >= min.ptos)
  {
    yy <- x[1:pos]
    xx <- r[1:pos]
    a1 <- cov(xx, yy) / var(xx)
    b1 <- mean(yy) - a1 * mean(xx)
  }

  a2 <- NULL
  b2 <- NULL

  if ((n - pos + 1) >= min.ptos)
  {
    yy <- x[pos:n]
    xx <- r[pos:n]
    a2 <- cov(xx,yy) / var(xx)
    b2 <- mean(yy) - a2 * mean(xx)
  }

  na <- as.double(NA)

  if (is.null(b1))
  {
    if (is.null(b2)) return(list(Tb = na, Tc = na, To = na, Ro = na,
                                 Tmin = min(x, na.rm = FALSE), Tmax = max(x, na.rm = FALSE),
                                 Nb = pos, Nc = n - pos + 1))
    return(list( Tb = na, Tc = b2, To = a2 * r[pos] + b2, Ro = r[pos],
                 Tmin=min(x, na.rm = FALSE), Tmax = max(x, na.rm = FALSE),
                 Nb = pos, Nc = n - pos + 1))
  } else
  {
    if (is.null(b2)) return(list(Tb = b1, Tc = na, To = a1 * r[pos] + b1,
                                 Ro = r[pos], Tmin = min(x, na.rm = FALSE),
                                 Tmax = max(x, na.rm = FALSE), Nb = pos, Nc = n - pos + 1))
    Ro = (b2 - b1) / (a1 - a2)
    return(list(Tb = b1, Tc = b2, To = a1 * Ro + b1, Ro = Ro,
                Tmin = min(x, na.rm = FALSE), Tmax = max(x, na.rm = FALSE),
                Nb = pos, Nc = n - pos + 1))
  }
}

maximizeR2 <- function(x, r, min.ptos = 3)
{
  n <- length(r)
  if (n < min.ptos) return(NA)
  opt <- as.integer(0)
  opt.val <- 0
  for (pos in min.ptos:(n - min.ptos + 1))
    tryCatch(
      {
        R2 <- 0

        yy <- x[1:pos]
        xx <- r[1:pos]
        cv <- cov(xx, yy)
        if (cv <= 0) next
        R2 <- (cv^2) / (var(xx) * var(yy))

        yy <- x[pos:n]
        xx <- r[pos:n]
        cv <- cov(xx, yy)
        if (cv >= 0) next
        R2 <- R2 + (cv^2) / (var(xx) * var(yy))

        if (R2 > opt.val)
        {
          opt <- pos
          opt.val <- R2
        }
      }, error = function(e) e)

  cv <- cov(x, r)
  R2 <- (cv^2) / (var(x) * var(r))
  if (cv > 0) if (R2 > opt.val) return(n)
  if (cv < 0) if (R2 > opt.val) return(1L)
  return(opt)
}

gaussfit <- function(ddd) tryCatch(
  {
    x <- ddd[["thetag"]]
    probit <- ddd[["probit"]]

    filter <- x > 0
    x <- x[filter]
    probit <- probit[filter]

    hlm <- lm(probit ~ x)
    b <- summary(hlm)$coefficients[1, 1]
    m <- summary(hlm)$coefficients[2, 1]

    theta50 <- -b/m
    sigma <- 1/m

    data.frame(theta50, sigma, R2 = summary(hlm)$r.squared)
  }, error = function(e) e)
