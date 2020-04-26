# Standardised measurement of seed functional traits: automated calculation of germination
# cardinal temperatures and thermal time using R 
# Authors; Emma Ladouceur, Hugh Pritchard, and Eduardo Fernandez-Pascual


# Beginners please refer to 'Supplementary File 3' for support in getting started

# Hash sign '#' is used for instructions and notes throughtout this script
# Lines that do not start with '#' are commands that you need to run
# To run a single line, left-click and position your I bar on the line, then hit run
# To run a part of a line or several lines, select and highlight them and hit run
# Running is done by shortkeys; Ctrl + Enter/R (RStudio/R, Windows) or Cmd + Return (Mac)
# Anytime you need help, type a '?' and the command next to it
# Try starting by running the command below with the shortcut key
?install.packages

# Install packages & dependencies (only the first time you use the script)
install.packages("binom", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("drc", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("plyr", dependencies = TRUE)
install.packages("segmented", dependencies = TRUE)

# Load the packages which you have just installed
# You can highlight them all and load them together, but check the output for errors
library(binom)
library(drc)
library(ggplot2)
library(plyr)
library(dplyr)
library(segmented)

# Set working directory
# Change this to fit your own working directory as detailed in 'Supplementary Information 3'
# Microsoft
setwd("C:/R")
# OSX
setwd("/Users/emmaladouceur/Desktop/Writing/8 Cardinal temps Seed Science research/")


# Load dataset, which is located in the working directory you just set
dat1 <- read.table("Supplementary information 1.txt", header = T)

# Exploring dataset; some basics
dat1                       # View the data in the output window
ncol(dat1)                 # How many columns does the dataset have?
nrow(dat1)                 # How many rows?
colnames(dat1)             # What are the column names?
head(dat1)                 # show the first 6 rows
head(dat1, n=20)           # show the first 20 rows
dat1[20:30,]               # show row number 20 to 30
dat1[,4:6]                 # show column number 4 to 6 
dat1[1,1]                  # show 1st column, 1st row item
dat1[,2]                   # What temperatures are there in the 2nd column 'Treatments'?
levels(dat1$Grouping)      # What are the 'levels' of the column, 'Grouping' 
summary(dat1)              # summarise all variables/columns
 
# IMPORTANT CHECKS
# Check that column names match the standard names, capitalization matters
# TRUE = OK; # FALSE = Revise column name
ifelse(colnames(dat1) == c("Grouping", "Treatment", "Dish", "Time", "G", "PG"), TRUE, FALSE)

# Treatment, Time, G and PG need to be numeric
# For example, letters in one cell of the Treatment column could render it a factor
# TRUE = OK; # FALSE = Revise dataset
is.numeric(dat1$Treatment)
is.numeric(dat1$Time)
is.numeric(dat1$G)
is.numeric(dat1$PG)
 
# Similarly, Grouping has to be a facto
# Using numbers for the level labels could make it numeric
# TRUE = OK; # FALSE = Revise dataset
is.factor(dat1$Grouping) 
   
# Ready to begin!


# STEP 1: Checking whether the data represents the full germination temperature range

# First filter data by the final scoring date
FSD <- ddply(dat1, .(Treatment, Dish, Grouping), filter, (Time == max(Time)))

# Write a function to estimate mean germination proportions 
# and the 95 % binomial confidence interval
# This is is a multi-line command, select all lines and run them
FGPfun <- function(x) {
          g <- sum(x$G)
          pg <- sum(x$PG)
          binom.confint(g, pg, methods = "logit")}

# Apply the function to every Treatment in each Grouping and store the output as 'FGP'
FGP <- ddply(FSD, .(Grouping, Treatment), FGPfun)

# You can view the output objects you have stored, such as 'FGP', by running their name
# These objects will last until you close the R session
FGP

# Export the new FGP dataset as a '.txt' file
# this will automatically be saved in the working directory you already set 
write.table(FGP, "Table S1 Final germination proportions.txt", sep = "\t", col.names = NA)

# Prepare to plot the final germination proportions
# Set the format for all future plots during this session
theme_update(panel.border = element_rect(linetype = "solid", colour = "black"))

# Create Figure 1 and call it 'FGPfig'
# You can make many aesthetic changes to these figures, see 'ggplot2' tutorials online
FGPfig <- ggplot(FGP, aes(y = mean, x = Treatment)) + 
          facet_wrap(~ Grouping) +
          geom_bar(stat = "identity", fill = "white", color = "black") +
          geom_errorbar(ymax = FGP$upper, ymin = FGP$lower, width = 1, lwd = 0.8) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                strip.background = element_rect(colour = "black", fill = "white")) +
          scale_y_continuous(labels = scales::percent) +
          coord_cartesian(ylim = c(0, 1)) +
          labs(x = "Temperature treatment (°C)", y = "Final germination")

# You could view Figure 1 by running 'FGPfig' (it would open in a new window)
# Now export Figure 1: 'FGPfig' as a '.tiff' file
# This will automatically be placed in your working directory
# Figures are saved with standard journal 2-column width, if you have many groupings they might be
# too small, and you may need to subset your dataset and make several figures
# You can adjust the width and height of the figures by changing the numbers after the arguments below
ggsave(filename = "Fig 1 Final germination percentage.tiff", plot = FGPfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


# STEP 2: Estimating germination rates from the cumulative germination curves

# Set a function to fit dose-response models to cumulative germination vs time
# These models will be used to estimate the time to succesive germination (e.g. t50)
# and the germination rates will be calculated as the inverse of time (e.g. r50 = 1/t50)
# Several dose-response models are available (see the documentation of the package 'drc')
# LL.2 is the log-logistic model (sometimes called Hill model)
# W1.2 and W2.2 are two types of Weibull models
# L.2 is a logistic model (sometimes called Boltzmann model)
# This function will select the best model for each treatment,
# based on Akaike's Information Criterion (AIC). Low AIC means the model fits the data better
GRfun <- function(x) {
         m1 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         m2 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = W1.2())
         m3 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = W2.2())
         m4 <- drm(G/PG ~ Time, data = x, type = "binomial", 
                   fct = L.3(fixed = c(NA, 1, NA)))
         e1 <- data.frame(ED(m1, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e2 <- data.frame(ED(m2, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e3 <- data.frame(ED(m3, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e4 <- data.frame(ED(m4, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F))
         e1$Model <- "LL.2"
         e2$Model <- "W1.2"
         e3$Model <- "W2.2"
         e4$Model <- "L.2"
         e1$AIC <- AIC(m1)
         e2$AIC <- AIC(m2)
         e3$AIC <- AIC(m3)
         e4$AIC <- AIC(m4)
         df1 <- rbind(e1, e2, e3, e4)
         colnames(df1)[1] <- "Times"
         df1$Rates <- 1/df1$Times
         df1$Deciles <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
         filter(df1, AIC == min(AIC))}

# Apply the above function to each combination of Grouping and Treatment
# Warnings may be printed when the function fails in some treatments 
# (e.g. because of low germination), ignore them
GR <- ddply(dat1, .(Grouping, Treatment), failwith(f = GRfun, quiet = T))
           
# View the table 'GR' with the germination rates 
GR

#For clarity, you can subset GR and view only the t50/r50
subset(GR, Deciles == "t50")

# Export Table 2 'GR' to a new file in your working directory
write.table(GR, "Table S2 Germination rates.txt", sep = "\t", col.names = NA)

# Now we'll plot the cumulative germinations and check function fit visually
# First you write a function to extract the parameters of the fitted functions
# Again, the function is selecting the dose-response model that fits better each case
CGfun <- function(x) {
         m1 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         m2 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = W1.2())
         m3 <- drm(G/PG ~ Time, data = x, type = "binomial", fct = W2.2())
         m4 <- drm(G/PG ~ Time, data = x, type = "binomial", 
                   fct = L.3(fixed = c(NA, 1, NA)))
         demo.fits <- expand.grid(conc = exp(seq(log(1.00e-04), 
                                  log(max(x$Time)), length = 100)))
         pm1 <- data.frame(predict(m1, newdata = demo.fits, interval = "confidence")) 
         pm2 <- data.frame(predict(m2, newdata = demo.fits, interval = "confidence")) 
         pm3 <- data.frame(predict(m3, newdata = demo.fits, interval = "confidence")) 
         pm4 <- data.frame(predict(m4, newdata = demo.fits, interval = "confidence")) 
         pm1$demo.fits.conc <- demo.fits$conc
         pm2$demo.fits.conc <- demo.fits$conc
         pm3$demo.fits.conc <- demo.fits$conc
         pm4$demo.fits.conc <- demo.fits$conc
         pm1$AIC <- AIC(m1)
         pm2$AIC <- AIC(m2)
         pm3$AIC <- AIC(m3)
         pm4$AIC <- AIC(m4)
         df2 <- rbind(pm1, pm2, pm3, pm4)
         filter(df2, AIC == min(AIC))}

# Then you apply the above function to each Grouping and Treatment
# The data point at time zero is modified slightly to improve plotting
# Warnings may be printed when the function fails in some treatments 
# (e.g. because of low germination), ignore them
CGfit <- ddply(dat1, .(Grouping, Treatment), failwith(f = CGfun, quiet = T))
dat1$Time1 <- dat1$Time
dat1$Time1[dat1$Time1 == 0] <- 1.00e-09

# Now you plot the figures according to the parameters you just extracted
# Again, figures might become too crowded if you have many groupings
# In the last line you'll see the label of the x axis (labs(x = "Time (days)")
# You may want to change the time units to the ones in your experiment
CGfig <- ggplot(CGfit, aes(x = demo.fits.conc, y = Prediction)) +
         geom_point(data = dat1, aes( x = Time1, y = G/PG), size = 0.8, alpha = 0.5)  +
         geom_line() +
         geom_ribbon(data = CGfit, aes(ymin = Lower, ymax = Upper), alpha = 0.2) + 
         facet_grid(Treatment ~ Grouping) +
         theme_bw() + 
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               strip.background = element_rect(colour = "black", fill = "white")) +
         scale_y_continuous(labels = scales::percent) +
         coord_cartesian(ylim = c(0, 1)) +
         labs(x = "Time (days)",y = "Final germination")

# You can view the figure by running 'CGfig', and check for weird fits and outliers
# Save Figure 2: CGfig to a '.tiff' file in your working directory
ggsave(filename = "Fig 2 Cumulative germination curves.tiff", plot = CGfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


# STEP 3a: Fitting a segmented model to the full germination temperature range 

# Set a function 'CTfun'  to fit segmented regression and calculate cardinal 
# temperatures and thermal time. The segmented function needs a preliminary
# estimate of the optimal temperature ('psi = to' in the 4th line), 
# so the function is giving it the temperature at which the germination rate is maximal
# Please note that each run of the segmented model can give slightly different fits
# and that some of the treatments may fail to fit in some of the runs
# Note also that the intercept of the second segment is calculated from the other
# model parameters and does not have a SE, t and p-value
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

# Apply the function 'CTfun' across 'Grouping' and 'Deciles' to calculate 
# cardinal temperatures and thermal time for each decil and call it 'CT'
# Warnings may be printed when the function fails in some deciles, ignore them
CT <- ddply(GR, .(Grouping, Deciles), failwith(f = CTfun, quiet = T))

# View Table 3: 'CT'
CT

#For clarity view only the calculations for the t50
subset(CT, Deciles == "t50")

# Write and export Table 3 to a file in your working directory
write.table(CT, "Table S3 Segmented model.txt", sep = "\t", col.names = NA)

# Now let's plot the segmented models
# Set a function to extract the parameters of the fits
BLfun <- function(x) {
         to <- filter(x, Rates == max(Rates))$Treatment
         lin.mod <- lm(Rates ~ Treatment, data = x)
         seg <- segmented(lin.mod, seg.Z = ~ Treatment, psi = to)
         data.frame(Treatment = unique(x)$Treatment, Rates = broken.line(seg)$fit)}
		 
# Apply the function to each Grouping and Decil
# Warnings may be printed when the function fails in some deciles, ignore them
BLfit <- ddply(GR, .(Grouping, Deciles), failwith(f = BLfun, quiet = T))

# Plot the germination rate vs germination temperature and the segmented models
# Again, in the last line you might want to change the time unit 
# of the y axis (days in this case) writing for example
# y = bquote("Germination rate ("*"hours"^-1 *")"))
BLfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point(size = 1) + 
         facet_grid(Grouping ~ Deciles, scales = "free_y") +
         geom_line(data = BLfit, colour = "black") +
         theme_bw() + 
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               strip.background = element_rect(colour="black", fill="white")) +
         labs(x = "Temperature treatment (°C)", 
		  y = bquote("Germination rate ("*"days"^-1 *")"))

# You can view the figure by running 'BLfig' to check for weird fits and outliers
# Outlying datapoints may suggest that some of the temperature treatments should be removed
# from the calculations
# Save Figure 3 to your working directory
ggsave(filename = "Fig 3 Germination rate vs T (segmented model).tiff", plot = BLfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


# STEP 3b: Fitting separate linear models to the sub- or 
# supraoptimal germination temperature range

# Set a function to fit 2 linear regressions and calculate cardinal T and thermal time
LMfun <- function(x) {
         tmax <- filter(x, Rates == max(Rates))$Treatment
         dsub <- filter(x, Treatment <= tmax)
         lin.sub <- lm(Rates ~ Treatment, data = dsub)
         a.sub <- lin.sub$coefficients[1]
         b.sub <- lin.sub$coefficients[2]        
         R.sub <- summary(lin.sub)$adj.r.squared
         dsupra <- filter(x, Treatment >= tmax)
         lin.supra <- lm(Rates ~ Treatment, data = dsupra)
         a.supra <- lin.supra $coefficients[1]
         b.supra <- lin.supra $coefficients[2]
         R.supra <- summary(lin.supra)$adj.r.squared
         data.frame(GRmax = filter(x, Rates == max(Rates))$Rates,
                    Treatment.max = tmax,
                    To = ifelse(b.sub > 0 & b.supra < 0,
                                (a.sub - a.supra) / (b.supra - b.sub), NaN),
                    n.sub = length(dsub$Treatment),
                    Intercept.sub = a.sub,
                    Intercept.sub.SE = summary(lin.sub)$coefficients[3],
                    Intercept.sub.t = summary(lin.sub)$coefficients[5],
                    Intercept.sub.p = summary(lin.sub)$coefficients[7],
                    Slope.sub = b.sub,
                    Slope.sub.SE = summary(lin.sub)$coefficients[4],
                    Slope.sub.t = summary(lin.sub)$coefficients[6],
                    Slope.sub.p = summary(lin.sub)$coefficients[8],
                    adj.R2.sub = R.sub,
                    Tb = ifelse(b.sub > 0, - a.sub / b.sub, NaN),
                    theta.sub = ifelse(b.sub > 0, abs(1/b.sub), NaN),
                    n.supra = length(dsupra$Treatment),
                    Intercept.supra = a.supra,
                    Intercept.supra.SE = summary(lin.supra)$coefficients[3],
                    Intercept.supra.t = summary(lin.supra)$coefficients[5],
                    Intercept.supra.p = summary(lin.supra)$coefficients[7],
                    Slope.supra = b.supra,
                    Slope.supra.SE = summary(lin.supra)$coefficients[4],
                    Slope.supra.t = summary(lin.supra)$coefficients[6],
                    Slope.supra.p = summary(lin.supra)$coefficients[8],
                    adj.R2.supra = R.supra,
                    Tc = ifelse(b.supra < 0, - a.supra / b.supra, NaN),
                    theta.supra = ifelse(b.supra < 0, abs(1/b.supra), NaN))}

# Apply the function to each Grouping and Decil
LM <- ddply(GR, .(Grouping, Deciles), failwith(f = LMfun, quiet = T))

# View the table 'LM'
LM

# View only the results for the t50
subset(LM, Deciles == "t50")

# Export Table 4 to your working directory
write.table(LM, "Table S4 Linear models.txt", sep = "\t", col.names = NA)

# Plot germination rate vs temperature and the linear models
# Again, consider changing the time unit of the y axis label in the last line
LM1 <- GR
LM2 <- LM
LM1$mb <- paste(LM1$Grouping, LM1$Deciles)
LM2$mb <- paste(LM2$Grouping, LM2$Deciles) 
LMfit <- merge(LM1, LM2, by = "mb")
LMfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         facet_grid(Grouping ~ Deciles, scales = "free_y") +         
         geom_point(size = 0.7)  +
         geom_abline(data = LM, aes(slope = Slope.sub, intercept = Intercept.sub)) +
         geom_abline(data = LM, aes(slope = Slope.supra, intercept = Intercept.supra)) +
         theme_bw() + 
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               strip.background = element_rect(colour = "black", fill = "white")) +
         labs(x = "Temperature treatment (°C)",
		  y = bquote("Germination rate ("*"days"^-1 *")"))

# You can view the figure by running 'LMfig' to check for weird fits and outliers
# Outlying datapoints may suggest that some of the temperature treatments should be removed
# from the calculations			  
# Export Table 4 to a file in your working directory
ggsave(filename = "Fig 4 Germination rate vs T (linear models).tiff", plot = LMfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


