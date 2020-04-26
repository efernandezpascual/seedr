# Load libraries

library(drc)
library(plyr)
library(segmented)
library(ggplot2)
library(dplyr)
library(binom)

# Load data

data <- read.table("data.txt", header = T)

### FINAL GERMINATION PROPORTIONS

# Filter data by final scoring date

FSD <- data %>% 
       group_by(Treatment, Dish) %>%
       filter(Time == max(Time))

# Write a function to estimate mean proportions and CI

FGPfun <- function(x) {
                       g <- sum(x$G)
                       pg <- sum(x$PG)
                       binom.confint(g, pg, 
                                     methods = "logit")
                       }

# Apply the function per treatment

FGP <- ddply(FSD, .(Treatment), FGPfun)

# Show and export final germination proportions

FGP
write.table(FGP, "c:/R/FGP.txt", sep = "\t", col.names = NA)

# Plot final germination proportions

FGPfig <- ggplot(FGP, aes(fill = Treatment, y = mean, x = Treatment)) +
          ylim(0.00, 1.00) + 
          geom_bar(stat = "identity") +
          geom_errorbar(ymax = FGP$upper, ymin = FGP$lower, width = 1, lwd = 1) +
          scale_fill_gradient(low = "blue", high = "red")
FGPfig

### GERMINATION TIME AND RATE

# Test the fit of different cumulative germination functions

FS <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                          LL2.5(), AR.2(), MM.2(), MM.3()), 
                          nested = FALSE, sorted = c("IC", "Res var", 
                          "Lack of fit", "no"), linreg = F, icfct = AIC)
      row.names(selected)[1]
                   }
ddply(data, .(Treatment), failwith(f = FS, quiet = T))

# Plot cumulative germination and check function fit visually

CGfun <- function(x) {
         fit <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         demo.fits <- expand.grid(conc = exp(seq(log(1.00e-04), 
                                  log(max(data$Time)), 
                                  length = 100))) 
         pm <- predict(fit, newdata = demo.fits, interval = "confidence")
         data.frame(demo.fits$conc, pm)    
                      }

CGfit <- ddply(data, .(Treatment), failwith(f = CGfun, quiet = T))

data$Time1 <- data$Time
data$Time1[data$Time1 == 0] <- 1.00e-09

CGfig <- ggplot(data, aes(x = Time1, y = G/PG)) +
         geom_point()  +
         geom_line(data = CGfit, 
                   aes( x = demo.fits.conc, y = Prediction)) +
         geom_ribbon(data = CGfit, 
                     aes(x = demo.fits.conc, y = Prediction, 
                     ymin = Lower, ymax = Upper), 
                     alpha = 0.2) + 
         facet_wrap(~ Treatment)

CGfig

# Set a function to fit the cumulative germination function 
# (set fct manually depending on previous step)

GRfun <- function(x) {
         fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         doses <- ED(fits, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F)
                      }

# Apply the function to each treatment and add decils and rates

GR <- ddply(data, .(Treatment), failwith(f = GRfun, quiet = T))
colnames(GR)[2] <- "Times"
GR$Rates <- 1/GR$Times
GR$Decils <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
                     
# Show and export germination times and rates

GR
write.table(GR, "c:/R/GR.txt", sep = "\t", col.names = NA)

### CARDINAL TEMPERATURES (IF BOTH SUB AND SUPRAOPTIMAL DATA AVAILABLE)

# Set a function to fit segmented regression and calculate cardinal temperatures

CTfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         segmented.mod <- segmented(lin.mod, seg.Z = ~ Treatment, psi = 20)
         a <- intercept(segmented.mod)
         asub <- a$Treatment[1, 1]
         asupra <- a$Treatment[2, 1]
         b <- slope(segmented.mod)
         bsub <- b$Treatment[1, 1]
         bsupra <- b$Treatment[2, 1]
         data.frame(Tb = - asub / bsub,
                    To = segmented.mod$psi[1, 2],
                    Tc = - asupra / bsupra)
                     }

# Calculate cardinal temperatures for each decil

CT <- ddply(GR, .(Decils), failwith(f = CTfun, quiet = T))

# Show and export cardinal temperatures

CT
write.table(CT, "c:/R/CT.txt", sep = "\t", col.names = NA)

### CARDINAL TEMPERATURES (IF ONLY SUBOPTIMAL DATA AVAILABLE)

# Set a function to fit segmented regression and calculate base temperature

TBfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         a <- lin.mod$coefficients[1]
         b <- lin.mod$coefficients[2]
         Tb <- (- a / b)
                     }

# Calculate cardinal temperatures for each decil

TB <- ddply(GR, .(Decils), failwith(f = TBfun, quiet = T))
colnames(TB)[2] <- "Tb"

# Show and export cardinal temperatures

TB
write.table(TB, "c:/R/TB.txt", sep = "\t", col.names = NA)

# Plot the GR versus time for each decil

TBfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point()  +
         geom_smooth(method = "lm", size = 1, se = T) +
         facet_wrap(~ Decils)
TBfig


### PART IN PROGRESS (SEGMENTED FIGURES)


lin.mod <- lm(t50~Treatment, data=times_wide)
segmented.mod <- segmented(lin.mod, seg.Z = ~Treatment, psi=20)
summary(segmented.mod)

#segmented
ggplot(GR, aes(x = Treatment, y = Rates)) +
geom_point() +
geom_line(data = segmented.mod) + 
facet_wrap(~ Decils)


qplot(Treatment, t50, group = Treatment > 27, geom = c('point', 'smooth'), 
      method = 'lm', se = F, data = times_wide)

