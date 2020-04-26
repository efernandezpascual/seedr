# Load libraries

library(drc)
library(plyr)
library(segmented)
library(ggplot2)
library(dplyr)
library(binom)

# Load data

data <- read.table("Carex diandra full range.txt", header = T)

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

FGP <- ddply(FSD, .(Grouping, Treatment), FGPfun)

# Show and export final germination proportions

FGP
write.table(FGP, "c:/R/Table 1 Final germination proportions.txt", sep = "\t", col.names = NA)

# Plot final germination proportions

FGPfig <- ggplot(FGP, aes(fill = Treatment, y = mean, x = Treatment)) + 
          facet_wrap(~ Grouping) +
          ylim(0.00, 1.00) + 
          geom_bar(stat = "identity") +
          geom_errorbar(ymax = FGP$upper, ymin = FGP$lower, width = 1, lwd = 1) +
          scale_fill_gradient(low = "blue", high = "red")
FGPfig
ggsave(filename = "Fig 1 Final germination proportions.tiff", plot = FGPfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)

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
ddply(data, .(Grouping, Treatment), failwith(f = FS, quiet = T))

# Plot cumulative germination and check function fit visually

CGfun <- function(x) {
         fit <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         demo.fits <- expand.grid(conc = exp(seq(log(1.00e-04), 
                                  log(max(data$Time)), 
                                  length = 100))) 
         pm <- predict(fit, newdata = demo.fits, interval = "confidence")
         data.frame(demo.fits$conc, pm)    
                      }

CGfit <- ddply(data, .(Grouping, Treatment), failwith(f = CGfun, quiet = T))

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
         facet_grid(Treatment ~ Grouping)

CGfig
ggsave(filename = "Fig 2 Cumulative germination curves.tiff", plot = CGfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)

# Set a function to fit the cumulative germination function 
# (set fct manually depending on previous step)

GRfun <- function(x) {
         fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         doses <- ED(fits, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F)
                      }

# Apply the function to each treatment and add decils and rates

GR <- ddply(data, .(Grouping, Treatment), failwith(f = GRfun, quiet = T))
colnames(GR)[3] <- "Times"
GR$Rates <- 1/GR$Times
GR$Decils <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
                     
# Show and export germination times and rates

GR
write.table(GR, "c:/R/Table 2 Germination rates.txt", sep = "\t", col.names = NA)

### CARDINAL TEMPERATURES (IF BOTH SUB AND SUPRAOPTIMAL DATA AVAILABLE)

# Set a function to fit segmented regression and calculate cardinal temperatures

CTfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         segmented.mod <- segmented(lin.mod, seg.Z = ~ Treatment, psi = 25)
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

CT <- ddply(GR, .(Grouping, Decils), failwith(f = CTfun, quiet = T))

# Show and export cardinal temperatures

CT
write.table(CT, "c:/R/Table 3 Cardinal temperatures.txt", sep = "\t", col.names = NA)

# Plot germination rates vs temperature per decil

BLfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         seg <- segmented(lin.mod, seg.Z = ~ Treatment, psi = 25)
         data.frame(Treatment = unique(x)$Treatment, 
                    Rates = broken.line(seg)$fit)
                     }

BLfit <- ddply(GR, .(Grouping, Decils), failwith(f = BLfun, quiet = T))

BLfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point() + 
         facet_grid(Grouping ~ Decils) +
         geom_line(data = BLfit, colour = "red")
BLfig

ggsave(filename = "Fig 3 Germination rate vs T (segmented model).tiff", plot = BLfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)

### CARDINAL TEMPERATURES (IF ONLY SUB OR SUPRAOPTIMAL DATA AVAILABLE)

# Set a function to fit segmented regression and calculate base or ceiling T

LMfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         a <- lin.mod$coefficients[1]
         b <- lin.mod$coefficients[2]
         - a / b
                     }

# Calculate cardinal temperatures for each decil

LM <- ddply(GR, .(Grouping, Decils), failwith(f = LMfun, quiet = T))
colnames(LM)[3] <- "Tb OR Tc"

# Show and export cardinal temperatures

LM
write.table(LM, "c:/R/Table 4 Base or ceiling temperatures.txt", sep = "\t", col.names = NA)

# Plot the GR versus time for each decil

LMfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point()  +
         geom_smooth(method = "lm", size = 1, se = T) +
         facet_grid(Grouping ~ Decils)
LMfig

ggsave(filename = "Fig 4 Germination rate vs T (linear model).tiff", plot = LMfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)

