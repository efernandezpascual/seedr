# Standardised measurement of seed functional traits: automated calculation of germination
# cardinal temperatures and thermal time using R 
# Authors; Emma Ladouceur, Hugh Pritchard, and Eduardo Fernandez-Pascual


# Beginners please refer to 'Supplementary File 3' for support in getting started

# Hash sign '#' is used for instructions and notes throughtout this script
# Lines that do not start with '#' are commands that you need to run
# To run a single line, left-click and position your I bar on the line, then hit run
# To run a part of a line or several lines, select and highlight them and hit run
# Running is done by shortkeys; Ctrl + Enter/Return (RStudio/R, Windows) or Cmd + Return (Mac)
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
setwd("/Users/emmaladouceur/Desktop/Academic/Manuscripts/2017_4 Cardinal temps/")


# Load dataset, which is located in the working directory you just set
dat <- read.table("S1 Example Data.txt", header = T)
# View the data in the output window by entering the next line
View(dat)

# Function to transform data into event-time format
TRAfun <- function(x) {
          data.frame(           
          Start = c(0, x$Time),
          End = c(x$Time, Inf),
          G = c(0, diff(x$G), max(x$PG) - tail(x$G, 1)),
          PG = max(x$PG)
          )}
dat1 <- ddply(subset(dat, Time != 0), .(Grouping, Treatment, Dish), TRAfun)
View(dat1)                   

#****************************
#STEPHANIE SUGGESTED THE SECTION BELOW
#WE NEED TO FILL IN # OF TREATMENTS AFTER ADDRESSING QUESTION
#ABOUT TEMPS MISSING FROM B (IF WE NEED TO?)
#************************************
# This dataset we have named 'dat1' and has 6 columns;
# Grouping: two 'factorial' 'levels': A & B which represent  two species 
# Temperature: ### treatments in intervals of 2.5 degrees Celcius
# Dish: represents number of  petri dishes, in this case there is one dish 
# Time; in days
# G: cumulative germination
# PG: number of viable seeds in each dish 

# Exploring dataset (dat1); some basics
ncol(dat1)                 # How many columns does the dataset have?
nrow(dat1)                 # How many rows?
colnames(dat1)             # What are the column names?
head(dat1)                 # show the first 6 rows as a default
head(dat1, n=20)           # show the first 20 rows
dat1[20:30,]               # show row number 20 to 30
dat1[,4:6]                 # show column number 4 to 6 
dat1[1,1]                  # show 1st column, 1st row item
levels(dat1$Grouping)      # What are the 'levels' of the column, 'Grouping' 
                           # The 'levels' command will only work for columns that are text,
                           # or 'factors' and will NOT work for numbers
summary(dat1)              # summarise all variables/columns
 
# IMPORTANT CHECKS
# Check that column names in the datafile match the names in the script, capitalization matters
# TRUE = OK; # FALSE = Revise column name
ifelse(colnames(dat1) == c("Grouping", "Treatment", "Dish", "Time", "G", "PG"), TRUE, FALSE)

# Treatment, Time, G and PG need to be numeric
# For example, letters in one cell of the Treatment column could render it a factor
# TRUE = OK; # FALSE = Revise dataset
is.numeric(dat1$Treatment)
is.numeric(dat1$Time)
is.numeric(dat1$G)
is.numeric(dat1$PG)

# Similarly, 'Grouping' must be text or a 'factor' and cannot be numeric
# TRUE = OK; # FALSE = Revise dataset
is.factor(dat1$Grouping) 
   
# Ready to begin!

# STEP 1: CHECK DATA FOR FULL TEMPERATURE RANGE

# First filter data by the final scoring date
FSD <- ddply(dat, .(Treatment, Dish, Grouping), filter, (Time == max(Time)))

# Write a function (FGPfun) to estimate mean final germination proportions  (FGP)
# and the 95 % binomial confidence interval
# This is a multi-line command, select all lines and run them together
FGPfun <- function(x) {
          g <- sum(x$G)
          pg <- sum(x$PG)
          binom.confint(g, pg, methods = "logit")}

# Apply the function to every Treatment in each Grouping and store the output as 'FGP'
FGP <- ddply(FSD, .(Grouping, Treatment), FGPfun)

# View the output objects you have stored, such as 'FGP', by running their name
# These objects will last until you close the R session
FGP

# Export the new FGP dataset as a '.txt' file
# this file will automatically be saved in the working directory you already set 
# but no output will appear in the output window
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
          labs(x = "Temperature treatment (?C)", y = "Final germination")

# You could view Figure 1 by running 'FGPfig' (it would open in a new window)
# Now export Figure 1: 'FGPfig' as a '.tiff' file
# This file will automatically be placed in your working directory
# Figures are saved with standard journal 2-column width, if you have many groupings the figures might be
# too small, and you may need to subset your dataset and make several figures
# If desired, adjust the width and height of the figures by changing the numbers after the arguments below
ggsave(filename = "Fig 1 Final germination percentage.tiff", plot = FGPfig,
       path = NULL, scale = 1, width = 70, height = 70, units = "mm", dpi = 300)


# STEP 2: ESTIMATE GERMINATION RATES FROM CURVES

# Write a function (GRfun) to fit time-to-event models to cumulative germination vs time
# These models will be used to estimate the time to germination (e.g. t50)
# and the germination rates will be calculated as the inverse of time (e.g. r50 = 1/t50)
# Two models are available and suitable for cumulative germination data (see the documentation of the package 'drc')
# The function below will select the best model for each treatment,
# based on Akaike's Information Criterion (AIC). Low AIC means the model fits the data better

# LL.2 is the log-logistic model (sometimes called Hill model)
# W2.2 is the Weibull model

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

# Now plot the cumulative germinations and check function fit visually
# First write a function to extract the parameters of the fitted cumulative germination functions
# Again, the function is selecting the time-to-event model that fits better each case
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

# Then apply the above function to each Grouping and Treatment
# The data point at time zero is modified slightly to improve plotting
# Warnings may be printed when the function fails in some treatments 
# (e.g. because of low germination), ignore them
CGfit <- ddply(dat1, .(Grouping, Treatment), failwith(f = CGfun, quiet = T))
dat$Time1 <- dat$Time
dat$Time1[dat$Time1 == 0] <- 1.00e-09

# Next plot the figures according to the parameters you just extracted
# Again, figures might become too crowded if you have many groupings
# In the last line you'll see the label of the x axis (labs(x = "Time (days)")
# Change the label  of for unit of time to what was used in your experiment
CGfig <- ggplot(CGfit, aes(x = demo.fits.conc, y = Prediction)) +
         geom_point(data = dat, aes( x = Time1, y = G/PG), size = 0.8, alpha = 0.5)  +
         geom_line() +
         geom_ribbon(data = CGfit, aes(ymin = Lower, ymax = Upper), alpha = 0.2) + 
         facet_grid(Grouping ~ Treatment) +
         theme_bw() + 
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               strip.background = element_rect(colour = "black", fill = "white")) +
         scale_y_continuous(labels = scales::percent) +
         coord_cartesian(ylim = c(0, 1)) +
         labs(x = "Time (days)",y = "Final germination")
CGfig
# You can view the figure by running 'CGfig', and check for weird fits and outliers
# Save Figure 2: CGfig to a '.tiff' file in your working directory
ggsave(filename = "Fig 2 Cumulative germination curves.tiff", plot = CGfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


# STEP 3A: CALCULATING CARDINAL TEMPERATURES (To,Tb, Tc)

#*******************************************
#** WHAT DOES PSI MEAN?
#*******************************************

# Set a function 'CTfun'  to fit segmented regression and calculate cardinal 
# temperatures and thermal time. The segmented function needs a preliminary
# estimate of the optimal temperature ('psi = to' in the 4th line), 
# so the function gives the temperature at which the germination rate is maximal
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

# Create Table 3 by applying the function 'CTfun' across 'Grouping' and 'Deciles' to calculate 
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
         new.x<- sort(c(unique(x)$Treatment, seg$psi[,2]))
         data.frame(Treatment=new.x,
                    Rates= predict(seg, newdata=data.frame(Treatment=new.x)))}
		 
# Apply the function (BLfit) to each Grouping and Decil
# Warnings may be printed when the function fails in some deciles, ignore them

#*******************************************
#**EXPLAIN WHT BL MEANS
#*******************************************

BLfit <- ddply(GR, .(Grouping, Deciles), failwith(f = BLfun, quiet = T))

# Plot the germination rate vs germination temperature and the segmented models
# Again, in the last line change the time unit  to the appropriate unit
# of the y axis (days in this case) writing for example
# y = bquote("Germination rate ("*"hours"^-1 *")"))
BLfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point(size = 1) + 
         facet_grid(Grouping ~ Deciles, scales = "free_y") +
         geom_line(data = BLfit, colour = "black") +
         theme_bw() + 
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               strip.background = element_rect(colour="black", fill="white")) +
         labs(x = "Temperature treatment (?C)", 
		  y = bquote("Germination rate ("*"days"^-1 *")"))
BLfig
# You can view the figure by running 'BLfig' to check for weird fits and outliers
# Outlying datapoints may suggest that some of the temperature treatments should be removed
# from the calculations
# Save Figure 3 to your working directory
ggsave(filename = "Fig 3 Germination rate vs T (segmented model).tiff", plot = BLfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


# STEP 3B: FIT SEPERARTE LINEAR MODELS TO SUB OR SUPRA RANGES

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
         labs(x = "Temperature treatment (?C)",
		  y = bquote("Germination rate ("*"days"^-1 *")"))

# You can view the figure by running 'LMfig' to check for weird fits and outliers
LMfig
# Outlying datapoints may suggest that some of the temperature treatments should be removed
# from the calculations			  
# Export Figure 4 to a file in your working directory
ggsave(filename = "Fig 4 Germination rate vs T (linear models).tiff", plot = LMfig,
       path = NULL, scale = 1, width = 173, height = 173, units = "mm", dpi = 300)


