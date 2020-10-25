#' Transforms dataset to physiodata format
#'
#' \code{physiodata} takes the user's dataset and transforms it to an object of
#' class "physiodata". This object will be used by the model-fitting functions,
#' and it can also be used to explore the data.
#'
#' @usage physiodata(d, t = "times", g = "germinated", pg = "germinable", x =
#'   "treatment", groups = NULL)
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
#' @return \code{physiodata} returns a S3 object of class "physiodata". The
#'   object is a list containing, for each group, treatment and scoring time:
#'   the cumulative germination count; the cumulative germination proportion;
#'   and the lower and upper bounds of the 95 % binomial confidence interval,
#'   calculated with the Wilson method as implemented in the package
#'   \code{binom}. The object can be used to explore the data using the generic
#'   functions \code{summary}, \code{barplot} and \code{plot}.
#' @examples
#' cent <- physiodata(centaury, x = "temperature")
#' cent
#' summary(cent) # average final germination proportions and germination rates per treatment
#' barplot(cent) # bar plots for the final germination proportions and germination rates
#' plot(cent) # cumulative germination curves
#' physiodata(grasses, x = "psi", groups = "species") # grouping dataset by species
#' @export
physiodata <- function(d, t = "times", g = "germinated", pg = "germinable",
                       x = "treatment", groups = NULL)
{
  dd <- data.table(d)
  setnames(dd, c(t, x), c("time", "treatment"))
  dd <- dd[, .(germinated = sum(get(g)), germinable = sum(get(pg))),
           by = c(groups, "treatment", "time")]
  dd[, germinable := max(germinable), by = c("treatment", groups)]
  setorderv(dd, c(groups, "treatment", "time"))
  dd[, cumulative := cumsum(get(g)), by = c("treatment", groups)]
  bci <- binom::binom.confint(dd$cumulative, dd$germinable, method = "wilson")
  dd <- cbind(dd, germination = bci[4:6])
  l <- list(proportions = dd, groups = groups)
  class(l) <- "physiodata"
  l
}

# physiodata generic functions

#' @export
print.physiodata <- function(x, ...)
{
  print(x$proportions)
}

#' @export
summary.physiodata <- function(object, ...)
{
  dd <- object$proportions[object$proportions[, .I[(time == max(time))], by = c(object$groups, "treatment")]$V1]
  dd[, c("germinated", "germinable", "cumulative") := NULL][]
  setorderv(dd, c(object$groups, "treatment"))
  dr <- object$proportions[, .(r50 = rates(d = .SD, fractions = 0.5, extrapolate.prange = 1)), by = c(object$groups, "treatment")]
  cbind(dd, dr[, list(r50)])
}

#' @export
barplot.physiodata <- function(height, ..., x.lab = "Treatment")
{
  dd <- summary(height)
  if(! is.null(height$groups)){
    listd <- split(dd, by = height$groups, drop = TRUE)
    ask.status <- par()$ask
    par(ask = TRUE)
    for(i in seq_along(listd)) {
      mfrow.status <- par()$mfrow
      oma.status <- par()$oma
      par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
      p <- barplot(listd[[i]]$germination.mean,
                   names.arg = as.numeric(listd[[i]][, treatment]),
                   ylim = c(0, 1),
                   xlab = x.lab,
                   ylab = "Final germination proportion")
      segments(p, listd[[i]]$germination.lower, p, listd[[i]]$germination.upper)
      arrows(p, listd[[i]]$germination.lower, p, listd[[i]]$germination.upper,
             lwd = 1.5, angle = 90, code = 3, length = 0.05)
      barplot(listd[[i]]$r50,
              names.arg = as.numeric(listd[[i]][, treatment]),
              xlab = x.lab,
              ylab = paste("Median germination rate"))
      mtext(names(listd)[i], line = 0, side = 3, outer = TRUE)
      par(mfrow = mfrow.status)
      par(oma = oma.status)
    }
    par(ask = ask.status)} else{
      mfrow.status <- par()$mfrow
      par(mfrow = c(1, 2))
      p <- barplot(dd$germination.mean,
                   names.arg = as.numeric(dd[, treatment]),
                   ylim = c(0, 1),
                   xlab = x.lab,
                   ylab = "Final germination proportion")
      segments(p, dd$germination.lower, p, dd$germination.upper)
      arrows(p, dd$germination.lower, p, dd$germination.upper,
             lwd = 1.5, angle = 90, code = 3, length = 0.05)
      barplot(dd$r50,
              names.arg = as.numeric(dd[, treatment]),
              xlab = x.lab,
              ylab = paste("Median germination rate"))
      par(mfrow = mfrow.status)
    }
}

#' @export
plot.physiodata <- function(x, ...)
{
  if(! is.null(x$groups)){
    listd <- split(x$proportions, by = x$groups, drop = TRUE)
    ask.status <- par()$ask
    par(ask = TRUE)
    for(i in seq_along(listd)) {
      colnumber <- listd[[i]][, .(n = length(unique(treatment)))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      X <- split(listd[[i]]$time, listd[[i]][, treatment])
      y <- split(listd[[i]]$germination.mean, listd[[i]][, treatment])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(X)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Germination proportion")
      mapply(lines, X, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(listd[[i]][, time]) + max(listd[[i]][, time])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(listd[[i]][, treatment], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)
      title(names(listd)[i])}
    par(ask = ask.status)} else{
      colnumber <- x$proportions[, .(n = length(unique(treatment)))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      X <- split(x$proportions$time, x$proportions[, treatment])
      y <- split(x$proportions$germination.mean, x$proportions[, treatment])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(X)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Germination proportion")
      mapply(lines, X, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(x$proportions[, time]) + max(x$proportions[, time])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(x$proportions[, treatment], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)}
}


# physiodata internal functions

rates <- function(d, fractions = (1:9)/10, extrapolate.prange = 1)
{
  pos <- c()
  pos0 <- c()
  for (i in 1:length(fractions))
  {
    posA <- match(FALSE, d$germination.mean < fractions[i], nomatch = NA)
    pos0A <- posA - 1
    if (is.na(posA))
    {
      posA <- length(d$germination.mean)
      pos0A <- match(FALSE, d$germination.mean < d$germination.mean[posA], nomatch = NA) - 1
    }
    if (pos0A == 0)
    {
      posA <- match(FALSE, d$germination.mean <= fractions[i], nomatch = NA)
      pos0A <- posA - 1
    }
    pos <- c(pos, posA)
    pos0 <- c(pos0, pos0A)
  }
  p <- d$germination.mean[pos]
  q <- d$germination.mean[pos0]
  y <- d$time[pos]
  x <- d$time[pos0]
  r <- x + (fractions - q) * (y - x) / (p - q)
  r[r > (extrapolate.prange * max(d$time))] <- NA
  1/r
}
