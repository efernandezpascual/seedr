setwd ("C:/EFP/Trabayu/Activos/2018 FICYT3/GT2 R tool/Bradford/Funciones Gil")
source("funciones.R")

bradtheta <- function(X, Y, Z) # Calculates Bradford's Theta
{
  Sxy <- cov(X, Y)
  Sxz <- cov(X, Z)
  Sxx <- var(X)
  Syz <- cov(Y, Z)
  Syy <- var(Y)
  Szz <- var(Z)

  a <- -2 * (Sxz ^ 2) * Syz + 2 * Sxy * Sxz * Szz
  b <- 2 * (Sxz ^ 2) * Syy + 4 * Sxy * Syz * Sxz - 2 * (Sxy ^ 2) * Szz - 4 * Syz * Sxy * Sxz
  c <- 2 * Syz * Sxy ^ 2 - 2 * Sxy * Sxz * Syy
 
  max(c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)))
}

braddf <- function(X) # Formats data frame
{ 
  X$cumulative = ave(X$germinated, X$temperature, X$psi, X$dish, FUN = cumsum) 
  X$germination = X$cumulative / X$germinable
  X$germination = ifelse(X$germination == 0, 0.001, X$germination) # No zeroes
  X$germination = ifelse(X$germination == 1, 0.999, X$germination) # No ones
  X$probit <- qnorm(X$germination, 0, 1)
  theta = bradtheta(X$probit, X$psi, 1/X$times)
  X$psibg = X$psi - theta / X$times
  data.frame(X)
}

bradtime <- function(X) # Calculates Bradford's hydrotime parameters
{
  bdf = braddf(X)
  blm = lm(probit ~ psibg, data = bdf)
  m = summary(blm)$coefficients[2, 1]
  b = summary(blm)$coefficients[1, 1]
  data.frame(theta = bradtheta(bdf$probit, bdf$psi, 1/bdf$times), 
             psib50 = -b/m, Sigma = 1/m, 
			 R2 = summary(blm)$r.squared)
}

bradplot <- function(X) # Plots Bradford's model
{
  bdf = braddf(X)
  plot(bdf$psibg, bdf$probit, xlab = "Psib(g)", ylab = "Probit germination")
  abline(lm(bdf$probit ~ bdf$psibg))
}

# One accession example

X <- read.table("Example.csv", header = T, sep = ",", dec = ".")

names(X) # Names need to comply!

names(X) <- c("temperature", "psi", "dish", 
              "times", "germinated", "germinable")

braddf(X)
bradtime(X)
bradplot(X)

# Multiple accessions example

X <- read.table("Hídrico.csv", header = T, sep = ",", dec = ".")

names(X) # Names need to comply!

names(X) <- c("species", "population", "temperature", "psi", "dish", 
              "times", "germinated", "germinable")
 
X <- split(X, X$species:X$population)

lapply(X, braddf)
lapply(X, bradtime)
par(mfrow = c(1, 2))
lapply(X, bradplot)

# ggplot to see differences amongst psis

X <- read.table("HÃ­drico.csv", header = T, sep = ",", dec = ".")

names(X) # Names need to comply!

names(X) <- c("species", "population", "temperature", "psi", "dish", 
              "times", "germinated", "germinable")

X %>% 
  #filter(population == "035AT") %>% 
  braddf %>%
  ggplot(aes(psibg, probit, color = as.factor(psi))) +
  geom_point() +
  geom_smooth(method = "lm")