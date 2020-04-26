

GRfun <- function(x) {
         m1 <- drm(G ~ Start + End, data = x, type = "event", fct = LL.2())
         m2 <- drm(G ~ Start + End, data = x, type = "event", fct = W2.2())
         e1 <- data.frame(ED(m1, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e2 <- data.frame(ED(m2, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e1$Model <- "LL.2"
         e2$Model <- "W2.2"
         e1$AIC <- AIC(m1)
         e2$AIC <- AIC(m2)
         df1 <- rbind(e1, e2)
         colnames(df1)[1] <- "Times"
         df1$Rates <- 1/df1$Times
         df1$Deciles <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
         filter(df1, AIC == min(AIC))}

CGfun <- function(x) {
         m1 <- drm(G ~ Start + End, data = x, type = "event", fct = LL.2())
         m2 <- drm(G ~ Start + End, data = x, type = "event", fct = W2.2())
         demo.fits <- expand.grid(conc = exp(seq(log(1.00e-04), 
                                  log(max(x$Start)), length = 100)))
         pm1 <- data.frame(predict(m1, newdata = demo.fits, interval = "confidence")) 
         pm2 <- data.frame(predict(m2, newdata = demo.fits, interval = "confidence")) 
         pm1$demo.fits.conc <- demo.fits$conc
         pm2$demo.fits.conc <- demo.fits$conc
         pm1$AIC <- AIC(m1)
         pm2$AIC <- AIC(m2)
         df2 <- rbind(pm1, pm2)
         filter(df2, AIC == min(AIC))}

CTfun <- function(x) {
         GRmax <- filter(x, Rates == max(Rates))$Rates
         tmax <- filter(x, Rates == max(Rates))$Treatment
         lin.mod <- lm(Rates ~ Treatment, data = x)
         segmented.mod <- segmented(lin.mod, seg.Z = ~ Treatment, psi = tmax)
         a <- intercept(segmented.mod)
         asub <- a$Treatment[1, 1]
         asupra <- a$Treatment[2, 1]
         b <- slope(segmented.mod)
         bsub <- b$Treatment[1, 1]
         bsupra <- b$Treatment[2, 1]
         data.frame(GRmax = GRmax,
                    Treatment.max = tmax,
                    To = segmented.mod$psi[1, 2],
                    n = length(x$Treatment),
                    adj.R2 = summary(segmented.mod)[9],
                    Intercept.sub = asub,
                    Intercept.sub.SE = summary(segmented.mod)$coefficients[5],
                    Intercept.sub.t = summary(segmented.mod)$coefficients[9],
                    Intercept.sub.p = summary(segmented.mod)$coefficients[13],
                    Slope.sub = bsub,
                    Slope.sub.SE = summary(segmented.mod)$coefficients[6],
                    Slope.sub.t = summary(segmented.mod)$coefficients[10],
                    Slope.sub.p = summary(segmented.mod)$coefficients[14],
                    Tb = ifelse(bsub > 0, - asub / bsub, NaN),
                    thetasub = ifelse(bsub > 0, 1/bsub, NaN),
                    Intercept.supra = asupra,
                    Slope.supra = bsupra,
                    Slope.supra.SE = slope(segmented.mod)$Treatment[4],
                    Slope.supra.t = summary(segmented.mod)$coefficients[11],
                    Slope.supra.p = summary(segmented.mod)$coefficients[15],
                    Tc = ifelse(bsupra < 0, - asupra / bsupra, NaN),
                    thetasupra = ifelse(bsupra < 0, -1/bsupra, NaN))}

BLfun <- function(x) {
         to <- filter(x, Rates == max(Rates))$Treatment
         lin.mod <- lm(Rates ~ Treatment, data = x)
         seg <- segmented(lin.mod, seg.Z = ~ Treatment, psi = to)
         new.x<- sort(c(unique(x)$Treatment, seg$psi[,2]))
         data.frame(Treatment=new.x,
                    Rates= predict(seg, newdata=data.frame(Treatment=new.x)))}

