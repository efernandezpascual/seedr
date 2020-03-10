#' @export

huidobro <- function(d, min.ptos = 3, method = c("Max R2","Max value"))
{
  dd <- d$rates[! is.na(d$rates$rate)]

  if ("Max value" %in% method)  # Tb<=Tmin, Tc>=Tmax for all quantiles
  {
    Ts <- dd[, cardinals(.SD[["treatment"]], rate, which.max(rate),
                          min.ptos = min.ptos), by = c(d$groups, "fraction")]
    Topt <- Ts[, .(Tb = min(c(mean(Tb, na.rm = TRUE), min(Tmin, na.rm = TRUE))),
                   Tc = max(c(mean(Tc, na.rm = TRUE), max(Tmax), na.rm = TRUE)),
                   To = mean(To, na.rm = TRUE)),
               by = c(d$groups)]
  }

  if ("Max R2" %in% method)  # Tb<=Tmin, Tc>=Tmax for all quantiles
  {
    Ts <- dd[, cardinals(.SD[["treatment"]], rate,
                          maximizeR2(.SD[["treatment"]], rate, min.ptos = min.ptos),
                          min.ptos = min.ptos), by = c(d$groups, "fraction")]
    Topt <- Ts[, .(Tb = min(c(mean(Tb, na.rm = TRUE), min(Tmin, na.rm = TRUE))),
                   Tc =max(c(mean(Tc, na.rm = TRUE), max(Tmax), na.rm = TRUE)),
                   To = mean(To, na.rm = TRUE)), by = c(d$groups)]
  }

  d$proportions <- Topt[d$proportions, on = c(d$groups)]  # Left join
  d$proportions[, suboptimal := (treatment <= To) * 1] # Suboptimal variable
  d$proportions[, thetag := time * (treatment - Tb) * suboptimal +
                  (Tc - treatment) * (1 - suboptimal), by = c(d$groups)] # thetag computation

  ddd <- split(d$proportions, d$proportions[, c(d$groups, "suboptimal"), with = FALSE])
  sapply(ddd, function(d) gaussfit(d))
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

gaussfit<-function(ddd)
{
  x <- ddd[["thetag"]]
  germination <- ddd[["germination.mean"]]

  filter <- x > 0
  x <- x[filter]
  germination <- germination[filter]

  probit = qnorm(germination, 0, 1)
  filter <- is.finite(probit)
  x <- x[filter]
  germination <- germination[filter]
  probit <- probit[filter]

  hlm <- lm.fit(matrix(c(x, rep(1, length(x))), ncol = 2), probit)
  m <- hlm$coefficients[1]
  b <- hlm$coefficients[2]

  theta50 <- -b/m
  sigma <- 1/m

  data.frame(theta50, sigma)
}
