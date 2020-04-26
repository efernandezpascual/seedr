library(drc)

data <- read.table("Carex diandra full range.txt", header = T)

fits <- drm(G/PG ~ Time, data = data, type = "binomial", fct = LL.2())

selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                    LL2.5(), AR.2(), MM.2(), MM.3()), 
                    nested = FALSE, sorted = c("IC", "Res var", 
                    "Lack of fit", "no"), linreg = F, icfct = AIC)

selected

row.names(selected)[1]
