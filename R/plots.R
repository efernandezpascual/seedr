# Plot germination proportions ==============================

#' @export

barplot.bradford <- function(d)
{
  p1 <- barplot(d$proportions$mean,
                names.arg = as.numeric(d$proportions[[1]]),
                ylim = c(0, 1),
                xlab = expression(paste(psi, " (MPa)")),
                ylab = "Final germination proportion")
  segments(p1, d$proportions$lower, p1, d$proportions$upper)
  arrows(p1, d$proportions$lower, p1, d$proportions$upper,
         lwd = 1.5, angle = 90, code = 3, length = 0.05)
}

#' @export

barplot.physiolist <- function(d)
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in d) {
    barplot(y)
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}

# Plot cumulative germination lines ==============================

#' @export

lines.bradford <- function(d)
{
  names(d$data)[1] <- "treatment"
  names(d$data)[2] <- "time"
  d1 <- d$data[, .(germinated = sum(cumulative), germinable = sum(germinable)),
               by = list(treatment, time)]
  d1[, proportion := germinated / germinable]
  setorder(d1, treatment, time)

  colnumber <- d1[, .(n = length(unique(treatment)))][[1]]
  colramp <- colorRampPalette(c("red", "orange", "yellow",
                                "green", "blue", "violet"))

  x <- split(d1$time, d1$treatment)
  y <- split(d1$proportion, d1$treatment)

  xpd.status <- par()$xpd
  mar.status <- par()$mar
  par(xpd = TRUE)
  par(mar = mar.status + c(0, 0, 0, 4))

  plot(1 : max(unlist(x)), ylim = (c(0, 1)), type = "n",
       xlab = "Time", ylab = "Final germination proportion")
  mapply(lines, x, y, col = colramp(colnumber), pch = 16, type = "o")
  legend(max(d1$time) + max(d1$time)*.05, .8,
         title = expression(paste(psi, " (MPa)")),
         legend = levels(as.factor(round(d$data[[1]], 1))), pch = 16,
         col = colramp(colnumber), lwd = 1, lty = 1)

  par(xpd = xpd.status)
  par(mar = mar.status)
}

#' @export

lines.physiolist <- function(d)
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in d) {
    lines(y)
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}

# Plot transformed data =======================================

#' @export

plot.bradford <- function(d, treatment) # Plots Bradford's model
{
  colnumber <- length(unique(d$data[[1]]))
  colramp <- colorRampPalette(c("red", "orange", "yellow",
                                "green", "blue", "violet"))
  d$data$color <- colramp(colnumber)[as.numeric(cut(d$data[[1]],breaks = colnumber))]

  plot(d$data$psibg, d$data$probit, col = d$data$color, pch = 16,
       xlab = expression(paste(psi[b], " (g)")), ylab = "Probit germination")
  abline(lm(d$data$probit ~ d$data$psibg))
  legend("topleft", title = expression(paste(psi, " (MPa)")),
         legend = levels(as.factor(round(d$data[[1]], 1))), pch = 16,
         col = colramp(colnumber))
}

#' @export

plot.physiolist <- function(d, treatment) # Plots Bradford's model
{
  ask.status <- par()$ask
  par(ask = TRUE)
  for(y in d) {
    plot(y)
    title(main = y$groups$group)
  }
  par(ask = ask.status)
}
