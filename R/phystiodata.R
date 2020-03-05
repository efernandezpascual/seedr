physiodata <- function(d, t = "times", g = "germinated", pg = "germinable", x = "temperature", reps = NULL, groups = NULL)
{
  dd <- data.table(d)
  dd <- dd[, .(germinated = sum(get(g)), germinable = sum(get(pg))),
           by = c(reps, groups, x, t)]
  dd[, germinable := max(germinable), by = c(x, reps, groups)]
  setorderv(dd, c(reps, groups, x, t))
  dd[, cumulative := cumsum(get(g)), by = c(x, reps, groups)]
  bci <- binom::binom.confint(dd$cumulative, dd$germinable, method = "wilson")
  dd <- cbind(dd, germination = bci[4:6])
  l <- list(data = dd, t = t, x = x, reps = reps, groups = groups)
  class(l) <- "physiodata"
  l
}

print.physiodata <- function(d)
{
  print(d$data)
}

summary.physiodata <- function(d)
{
  dd <- d$data[(get(d$t) == max(get(d$t))),
               list(times = get(d$t), germination.mean, germination.lower, germination.upper),
               by = c(d$groups, d$x, d$reps)]

  dd <- d$data[d$data[, .I[(get(d$t) == max(get(d$t)))], by = c(d$groups, d$x, d$reps)]$V1]

  dd <- dd[, list(species = get(d$group), treatment = get(d$x), reps = get(d$reps), times = get(d$t), germination.mean, germination.lower, germination.upper)]


  setorderv(dd, c(d$groups, d$x, d$reps))
  dd
}

barplot.physiodata <- function(d, x.lab = "Treatment")
{
  dd <- summary(d)
  if(! is.null(d$groups)){
  listd <- split(dd, by = d$groups, drop = TRUE)
  ask.status <- par()$ask
  par(ask = TRUE)
  for(i in seq_along(listd)) {
    p <- barplot(listd[[i]]$germination.mean,
                 names.arg = as.numeric(listd[[i]][, get(d$x)]),
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
                 names.arg = as.numeric(dd[, get(d$x)]),
                 ylim = c(0, 1),
                 xlab = x.lab,
                 ylab = "Final germination proportion")
    segments(p, dd$germination.lower, p, dd$germination.upper)
    arrows(p, dd$germination.lower, p, dd$germination.upper,
           lwd = 1.5, angle = 90, code = 3, length = 0.05)
  }
}

plot.physiodata <- function(d)
{
  if(! is.null(d$groups)){
    listd <- split(d$data, by = d$groups, drop = TRUE)
    ask.status <- par()$ask
    par(ask = TRUE)
    for(i in seq_along(listd)) {
      colnumber <-listd[[i]][, .(n = length(unique( get(d$x))))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      x <- split(listd[[i]]$time, listd[[i]][, get(d$x)])
      y <- split(listd[[i]]$germination.mean, listd[[i]][, get(d$x)])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(x)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Final germination proportion")
      mapply(lines, x, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(listd[[i]][, get(d$t)]) + max(listd[[i]][, get(d$t)])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(listd[[i]][, get(d$x)], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)
      title(names(listd)[i])}
    par(ask = ask.status)} else{
      colnumber <- d$data[, .(n = length(unique( get(d$x))))][[1]]
      colramp <- colorRampPalette(c("violet", "blue", "green",
                                    "yellow", "orange", "red"))
      x <- split(d$data$time, d$data[, get(d$x)])
      y <- split(d$data$germination.mean, d$data[, get(d$x)])
      xpd.status <- par()$xpd
      mar.status <- par()$mar
      par(xpd = TRUE)
      par(mar = mar.status + c(0, 0, 0, 4))
      plot(1 : max(unlist(x)), ylim = (c(0, 1)), type = "n",
           xlab = "Time", ylab = "Final germination proportion")
      mapply(lines, x, y, col = colramp(colnumber), pch = 16, type = "o")
      legend(max(d$data[, get(d$t)]) + max(d$data[, get(d$t)])*.05, 1.1,
             title = "Treatment",
             legend = levels(as.factor(round(d$data[, get(d$x)], 1))), pch = 16,
             col = colramp(colnumber), lwd = 1, lty = 1)
      par(xpd = xpd.status)
      par(mar = mar.status)}
}

