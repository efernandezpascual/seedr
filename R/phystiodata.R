#' @export

physiodata <- function(d, t = "times", g = "germinated", pg = "germinable", x = "temperature", groups = NULL)
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
  l <- list(data = dd, groups = groups)
  class(l) <- "physiodata"
  l
}

# physiodata generic functions

#' @export

print.physiodata <- function(d)
{
  print(d$data)
}

#' @export

summary.physiodata <- function(d)
{
  dd <- d$data[d$data[, .I[(time == max(time))], by = c(d$groups, "treatment")]$V1]
  dd[, c("germinated", "germinable", "cumulative") := NULL][]
  setorderv(dd, c(d$groups, "treatment"))
  dd
}

#' @export

barplot.physiodata <- function(d, x.lab = "Treatment")
{
  dd <- summary(d)
  if(! is.null(d$groups)){
  listd <- split(dd, by = d$groups, drop = TRUE)
  ask.status <- par()$ask
  par(ask = TRUE)
  for(i in seq_along(listd)) {
    p <- barplot(listd[[i]]$germination.mean,
                 names.arg = as.numeric(listd[[i]][, treatment]),
                 ylim = c(0, 1),
                 xlab = x.lab,
                 ylab = "Final germination proportion")
    segments(p, listd[[i]]$germination.lower, p, listd[[i]]$germination.upper)
    arrows(p, listd[[i]]$germination.lower, p, listd[[i]]$germination.upper,
           lwd = 1.5, angle = 90, code = 3, length = 0.05)
    title(names(listd)[i])
    }
  par(ask = ask.status)} else{
    p <- barplot(dd$germination.mean,
                 names.arg = as.numeric(dd[, treatment]),
                 ylim = c(0, 1),
                 xlab = x.lab,
                 ylab = "Final germination proportion")
    segments(p, dd$germination.lower, p, dd$germination.upper)
    arrows(p, dd$germination.lower, p, dd$germination.upper,
           lwd = 1.5, angle = 90, code = 3, length = 0.05)
  }
}

#' @export

plot.physiodata <- function(d)
{
  if(! is.null(d$groups)){
    listd <- split(d$data, by = d$groups, drop = TRUE)
    ask.status <- par()$ask
    par(ask = TRUE)
    for(i in seq_along(listd)) {
      colnumber <- listd[[i]][, .(n = length(unique(treatment)))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      x <- split(listd[[i]]$time, listd[[i]][, treatment])
      y <- split(listd[[i]]$germination.mean, listd[[i]][, treatment])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(x)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Final germination proportion")
      mapply(lines, x, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(listd[[i]][, time]) + max(listd[[i]][, time])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(listd[[i]][, treatment], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)
      title(names(listd)[i])}
    par(ask = ask.status)} else{
      colnumber <- d$data[, .(n = length(unique(treatment)))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      x <- split(d$data$time, d$data[, treatment])
      y <- split(d$data$germination.mean, d$data[, treatment])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(x)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Final germination proportion")
      mapply(lines, x, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(d$data[, time]) + max(d$data[, time])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(d$data[, treatment], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)}
}

