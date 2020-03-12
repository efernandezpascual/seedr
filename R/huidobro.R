#' @export

huidobro <- function(d, min.ptos = 3, tops = c("Max R2","Max value"),
                     fractions = (1:9)/10, extrapolate.prange = 1)
{
  rates <- d[, .(fraction = fractions, rate = rates(d = .SD, fractions = fractions, extrapolate.prange = extrapolate.prange)), by = c("treatment")]
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
}

#' @export

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

#' @export

maximizeR2 <- function(x, r, min.ptos = 3)
{
  n <- length(r)
  if (n < min.ptos) return(NA)
  opt <- as.integer(0)
  opt.val <- 0
  for (pos in min.ptos:(n - min.ptos + 1))
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
  }

  cv <- cov(x, r)
  R2 <- (cv^2) / (var(x) * var(r))
  if (cv > 0) if (R2 > opt.val) return(n)
  if (cv < 0) if (R2 > opt.val) return(1L)
  return(opt)
}

#' @export

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

#' @export

print.huidobro <- function(y)
{
  if(is.nan(y$parameters$cardinals$Tc)) {
    cat("Garcia-Huidobro's suboptimal thermal time model", "\n",
        "Suboptimal temperature levels in experiment:",
        round(y$parameters$levels$suboptimal, 1), "\n",
        "Tb - Base temperature:",
        round(y$parameters$cardinals$Tb, 1), "\n",
        "theta50 1 - Suboptimal thermal time (median):",
        round((y$parameters$suboptimal$theta50), 2), "\n",
        "Sigma of the suboptimal thermal time:",
        round((y$parameters$suboptimal$sigma), 2), "\n",
        "R2 of the suboptimal model:",
        round(y$parameters$suboptimal$R2, 2), "\n")} else {
          if(is.nan(y$parameters$cardinals$To)) {
            cat("Garcia-Huidobro's supraoptimal thermal time model", "\n",
                "Supraoptimal temperature levels in experiment:",
                round(y$parameters$levels$supraoptimal, 1), "\n",
                "Tc - Ceiling temperature:",
                round(y$parameters$cardinals$Tc, 1), "\n",
                "theta50 1 - Suboptimal thermal time (median):",
                round(y$parameters$supraoptimal$theta50, 2), "\n",
                "Sigma of the suboptimal thermal time:",
                round(y$parameters$supraoptimal$sigma, 2), "\n",
                "R2 of the suboptimal model:",
                round(y$parameters$supraoptimal$R2, 2), "\n", "\n")} else
                  cat("Garcia-Huidobro's thermal time model", "\n",
                      "Suboptimal temperature levels in experiment:",
                      round(y$parameters$levels$suboptimal, 1), "\n",
                      "Supraoptimal temperature levels in experiment:",
                      round(y$parameters$levels$supraoptimal, 1), "\n",
                      "Tb - Base temperature:",
                      round(y$parameters$cardinals$Tb, 1), "\n",
                      "To - Optimal temperature:",
                      round(y$parameters$cardinals$To, 1), "\n",
                      "Tc - Ceiling temperature:",
                      round(y$parameters$cardinals$Tc, 1), "\n",
                      "theta50 1 - Suboptimal thermal time (median):",
                      round((y$parameters$suboptimal$theta50), 2), "\n",
                      "Sigma of the suboptimal thermal time:",
                      round((y$parameters$suboptimal$sigma), 2), "\n",
                      "R2 of the suboptimal model:",
                      round(y$parameters$suboptimal$R2, 2), "\n",
                      "theta50 2 - Supraoptimal thermal time (median):",
                      round(y$parameters$supraoptimal$theta50, 2), "\n",
                      "Sigma of the supraoptimal thermal time:",
                      round(y$parameters$supraoptimal$sigma, 2), "\n",
                      "R2 of the supraoptimal model:",
                      round(y$parameters$supraoptimal$R2, 2), "\n", "\n")}
}

#' @export

summary.huidobro <- function(y)
{
  data.table(
    method = class(y),
    n.treatments.sub = length(y$parameters$levels$suboptimal),
    n.treatments.sup = length(y$parameters$levels$supraoptimal),
    Tb = y$parameters$cardinals$Tb,
    To = y$parameters$cardinals$To,
    Tc = y$parameters$cardinals$Tc,
    theta50.sub = y$parameters$suboptimal$theta50,
    sigma.sub = y$parameters$suboptimal$sigma,
    R2.sub = y$parameters$suboptimal$R2,
    theta50.sup = y$parameters$supraoptimal$theta50,
    sigma.sup = y$parameters$supraoptimal$sigma,
    R2.sup = y$parameters$supraoptimal$R2)
}

#' @export

plot.huidobro <- function(d) # Plots Bradford's model
{
  if(is.nan(d$parameters$cardinals$Tc)) {
  xpd.status <- par()$xpd
  mar.status <- par()$mar

  par(mar = mar.status + c(0, 0, 0, 4))
  colnumber <- length(unique(d$data$suboptimal$treatment))
  colramp <- colorRampPalette(c("red", "orange", "yellow",
                                "green", "blue", "violet"))
  d$data$suboptimal$color <- colramp(colnumber)[as.numeric(cut(d$data$suboptimal$treatment, breaks = colnumber))]

  plot(d$data$suboptimal$thetag, d$data$suboptimal$probit, col = d$data$suboptimal$color, pch = 16,
       xlab = expression(paste("Suboptimal ", theta, " (g)")), ylab = "Probit germination")
  abline(lm(d$data$suboptimal$probit ~ d$data$suboptimal$thetag))
  par(xpd = TRUE)
  legend(max(d$data$suboptimal[, thetag]) + max(d$data$suboptimal[, thetag])*.05, 2,
         title = "ºC",
         legend = levels(as.factor(round(d$data$suboptimal$treatment, 1))), pch = 16,
         col = colramp(colnumber))
  par(xpd = xpd.status)
  par(mar = mar.status)} else
  {if(is.nan(d$parameters$cardinals$To)) {
    xpd.status <- par()$xpd
    mar.status <- par()$mar

    par(mar = mar.status + c(0, 0, 0, 4))
    colnumber <- length(unique(d$data$supraoptimal$treatment))
    colramp <- colorRampPalette(c("red", "orange", "yellow",
                                  "green", "blue", "violet"))
    d$data$supraoptimal$color <- colramp(colnumber)[as.numeric(cut(d$data$supraoptimal$treatment, breaks = colnumber))]

    plot(d$data$supraoptimal$thetag, d$data$supraoptimal$probit, col = d$data$supraoptimal$color, pch = 16,
         xlab = expression(paste("Supraoptimal ", theta, " (g)")), ylab = "Probit germination")
    abline(lm(d$data$supraoptimal$probit ~ d$data$supraoptimal$thetag))
    par(xpd = TRUE)
    legend(max(d$data$supraoptimal[, thetag]) + max(d$data$supraoptimal[, thetag])*.05, 2,
           title = "ºC",
           legend = levels(as.factor(round(d$data$supraoptimal$treatment, 1))), pch = 16,
           col = colramp(colnumber))
    par(xpd = xpd.status)
    par(mar = mar.status)} else {
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      ask.status <- par()$ask
      par(ask = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))

      colnumber <- length(unique(d$data$suboptimal$treatment))
      colramp <- colorRampPalette(c("red", "orange", "yellow",
                                    "green", "blue", "violet"))
      d$data$suboptimal$color <- colramp(colnumber)[as.numeric(cut(d$data$suboptimal$treatment, breaks = colnumber))]

      plot(d$data$suboptimal$thetag, d$data$suboptimal$probit, col = d$data$suboptimal$color, pch = 16,
           xlab = expression(paste("Suboptimal ", theta, " (g)")), ylab = "Probit germination")
      abline(lm(d$data$suboptimal$probit ~ d$data$suboptimal$thetag))
      par(xpd = TRUE)
      legend(max(d$data$suboptimal[, thetag]) + max(d$data$suboptimal[, thetag])*.05, 2,
             title = "ºC",
             legend = levels(as.factor(round(d$data$suboptimal$treatment, 1))), pch = 16,
             col = colramp(colnumber))
      par(xpd = xpd.status)

      colnumber <- length(unique(d$data$supraoptimal$treatment))
      colramp <- colorRampPalette(c("red", "orange", "yellow",
                                    "green", "blue", "violet"))
      d$data$supraoptimal$color <- colramp(colnumber)[as.numeric(cut(d$data$supraoptimal$treatment, breaks = colnumber))]

      plot(d$data$supraoptimal$thetag, d$data$supraoptimal$probit, col = d$data$supraoptimal$color, pch = 16,
           xlab = expression(paste("Supraoptimal ", theta, " (g)")), ylab = "Probit germination")
      abline(lm(d$data$supraoptimal$probit ~ d$data$supraoptimal$thetag))
      par(xpd = TRUE)
      legend(max(d$data$supraoptimal[, thetag]) + max(d$data$supraoptimal[, thetag])*.05, 2,
             title = "ºC",
             legend = levels(as.factor(round(d$data$supraoptimal$treatment, 1))), pch = 16,
             col = colramp(colnumber))

      par(xpd = xpd.status)
      par(mar = mar.status)
      par(ask = ask.status)
    }}
}

