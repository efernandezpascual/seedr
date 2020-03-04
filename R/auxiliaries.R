# Format input data.frame =========================================

#' @export

physiotable <- function(d, t, g, pg, x, reps = NULL, groups = NULL)
{
  d <- data.table(d)
  d <- d[, .(germinated = sum(get(g)), germinable = sum(get(pg))),
         by = c(x, t, groups, reps)]
  d[, germinable := max(germinable), by = c(x, reps, groups)]
  setorderv(d, c(groups, reps, x, t))
  d[, cumulative := cumsum(germinated), by = c(x, reps, groups)]
  bci <- binom.confint(d$cumulative, d$germinable, method = "wilson")
  cbind(d, bci[4:6])
}

# Calculate germination proportions with CI ==========================

#' @export
#' @import binom

proportions <- function(d, x, reps = NULL) # This function calculates binomials!
{
  d1 <- d[, .(germinated = max(cumulative), germinable = max(germinable)), by = c(x, reps)]
  d1 <- d1[, .(germinated = sum(germinated), germinable = sum(germinable)), by = x]
  bci <- binom.confint(d1$germinated, d1$germinable, method = "wilson")
  cbind(d1, bci[4:6])
}
