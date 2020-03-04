# Print generics =========================================

#' @export

print.bradford <- function(y)
{
  cat("Bradford's hydrotime model", "\n",
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

#' @export

print.physiolist <- function(d)
{
  cat("A list of physiological time germination models", "\n",
      "calculated for the following", length(d), "groups:", "\n", "\n")
  for(y in d) {
    cat(y$groups$group, "\n")
    print(y)
  }
}

# Summary generics show germination proportions ==========================

#' @export

summary.bradford <- function(y) y$proportions

#' @export

summary.physiolist <- function(d) lapply(d, function(y) y$proportions)

# Coefficient generics show model parameters =============================

#' @export

coef.bradford <- function(y)
{
  data.table(
    method = class(y),
    n.viable = y$parameters$n.viable,
    n.germinated = y$parameters$n.germinated,
    n.treatments = length(y$parameters$levels),
    theta = y$parameters$theta,
    psib50 = y$parameters$psib50,
    sigma = y$parameters$sigma,
    R2 = y$parameters$R2)
}

#' @export

coef.physiolist <- function(d)
{
  l <- lapply(d, function(y) data.frame(
    y$groups[,-ncol(y$groups), with = FALSE], coef(y)))
  rbindlist(l)
}
