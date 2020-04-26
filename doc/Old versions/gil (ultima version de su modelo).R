

#' Get sample distribution functions for each psi based on sampled data
#' @import data.table
# @export

GerminationTT2 <- function(ddt)
{
  x <- as.data.table(ddt)
  x <- x[, sumacum := cumsum(germinated), by = list(dish, psi)]
  y <- x[, .(germinable = max(germinable)), by = list(dish, psi)]
  germinable <- y[, .(germinable = sum(germinable)), by = psi]
  x <- x[, rbind(.(times = c(0, 0), sumacum = c(0, 0),
                   germinable = c(NA, NA)), .SD), by = list(dish, psi)] #Tuve que añadir fill = T
  x <- x[, .(times = tail(times, .N - 1), ni = tail(sumacum, .N -1) -
               head(sumacum, .N - 1)), by = list(dish, psi)]
  x <- x[, .(ni = sum(ni), Prob = germinable$germinable[germinable$psi == psi]),
         by = list(times, psi)]
  x[, .(times, ni, Prob = cumsum(ni) / Prob), by = list(psi)]
}

# @export
mapprox <- function(x, y, xout)
{
  if (max(x) - min(x) < 1E-10)
    return(list(x = xout, y = rep(as.numeric(NA), times = length(xout))))
  approx(x, y, xout)
}

#' Estimations of Theta
# @export
ModelBasedThetaV4 <- function(it)
{
  #   it<-it[it$Prob>0,]
  #   it<-it[it$ni>0,]
  ii <- it[, cbind(mapprox(.SD$Prob, .SD$times, it$Prob),
                   data.frame(times = it$times, Treat = it$psi, ni = it$ni)),
           by = psi]
  ii <- ii[ii$psi > ii$Treat, ]
  ratio <- (1 - ii$y/ii$times)
  ii$Theta <- (ii$psi - ii$Treat) * ii$y / ratio
  #   ii[is.finite(ii$Theta)&(ii$Theta>0)&(ii$ni>0)&(ratio>1E-10)]
  ii[is.finite(ii$Theta) & (ii$ni > 0) & (ratio > 1E-10)]
}
