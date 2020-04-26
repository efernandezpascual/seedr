# Load libraries
library(drc)
library(plyr)
library(segmented)

# Load data
data <- read.table("data.txt", header = T)

# Test the fit of different cumulative germination functions
fs <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                          LL2.5(), AR.2(), MM.2(), MM.3()), 
                          nested = FALSE, sorted = c("IC", "Res var", 
                          "Lack of fit", "no"), linreg = F, icfct = AIC)
      row.names(selected)[1]
      }
ddply(data, .(Treatment), failwith(f = fs, quiet = T))

# Set a function to fit the cumulative germination function 
# (set fct manually depending on previous step)
tg <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      doses <- ED(fits, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F)
      }

# Apply the function to each treatment and add deciland rates
times <- ddply(data, .(Treatment), failwith(f = tg, quiet = T))
times$decils <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
times$rates <- 1/times$Estimate

# Set a function to fit segmented regression and calculate cardinal temperatures
seg <- function(x) {
       lin.mod <- lm(rates ~ Treatment, data = x)
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
Cardinal.T <- ddply(times, .(decils), failwith(f = seg, quiet = T))
Cardinal.T