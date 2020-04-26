
# An easy and automated calculation of the germination cardinal temperatures and thermal time using R
# Authors; Emma Ladouceur, Hugh Pritchard, and Eduardo Fernandez-Pascual

# Beginners please refer to 'Supplementary File 3' for support in getting started

# Hash sign '#' is used for instructions and notes throughtout this script
# Shortcut keys are recommended to be used; Ctrl + Enter/R (RStudio/R, Windows) or Cmd + Return (Mac)
# You don't have to highlight the line to run it, simply left-click and position your I bar on the line, then hit run
# Anytime you need help, type a '?' and the command next to it
# Try starting by running the command below with the shortcut key
?install.packages

# Install packages & dependencies
install.packages("binom", dependencies=TRUE)
install.packages("dplyr", dependencies=TRUE)
install.packages("drc", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("plyr", dependencies=TRUE)
install.packages("segmented", dependencies=TRUE)

# Load required libraries which you have just installed
library(binom)
library(dplyr)
library(drc)
library(ggplot2)
library(plyr)
library(segmented)

# Set working directory
# Change this to fit your own working directory as detailed in 'Supplementary Information 3'
setwd("/Users/emmaladouceur/Desktop/Writing/8 Cardinal temps Seed Science research")

# Load dataset, which is located in the working directory you just set
data <- read.table("Supplementary Information  1_1_17_EL.txt", header = T)

# Exploring data; some basics
View(data)                 # View the data in a seperate window
data                       # or directly in the output window
ncol(data)                 # How many columns does the dataset have?
nrow(data)                 # How many rows?
colnames(data)             # What are the column names?
head(data)                 # show the first 6 rows
head(data, n=20)           # show the first 20 rows
data[20:30,]               # show row number 20 to 30
data[,4:6]                 # show column number 4 to 6 
data[1,1]                  # show 1st column, 1st row item
data[,2]                   # What temperatures are there in the 2nd column 'Treatments'?
is.numeric(data$Treatment) # is 'Treatment' data numeric? TRUE=YES, FALSE=NO
is.factor(data$Grouping)   # is 'Grouping' data factorial?
levels(data$Grouping)      # What are the 'levels' of the column, 'Grouping' 
summary(data)              # summarise all variables/columns
    
# Ready to begin!

# Step 1: Checking whether the data represents the full germination temperature range

# First filter data by the final scoring date
# this is is a multi-line command, place your cursor at the end of it and use your shortcut key to enter it
FSD <- data %>% 
       group_by(Treatment, Dish) %>%
       filter(Time == max(Time))

# A function to estimate mean germination proportions and the associated confidence interval
FGPfun <- function(x) {
                       g <- sum(x$G)
                       pg <- sum(x$PG)
                       binom.confint(g, pg, 
                                     methods = "logit")
                       }

# Apply the function per for every treatment in each grouping and call it 'FGP'
FGP <- ddply(FSD, .(Grouping, Treatment), FGPfun)

# View your first table, Table 1: FGP
View(FGP)
# write the new dataset Table 1: 'FGP' to a '.txt' file, seperate the fields (sep), with a tab ("\t")
# this will automatically be written to the working directory you already set 
write.table(FGP, "Table 1 Final germination proportions.txt", sep = "\t", col.names = NA)

# Prepare to plot final germination proportions
# Set the format for all future plots, nothing will happen here, its purely settings for future action
theme_update(panel.border = element_rect(linetype = "solid", colour = "black"))

# Create Figure 1 and call it 'FGPfig'
FGPfig <- ggplot(FGP, aes(y = mean, x = Treatment)) + 
          facet_wrap(~ Grouping) +
          ylim(0.00, 1.00) + 
          geom_bar(stat = "identity", fill= "white", color="black") +
          geom_errorbar(ymax = FGP$upper, ymin = FGP$lower, width = 1, lwd = 0.8) + 
          theme_bw() + theme(panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             strip.background = element_rect(colour="black", fill="white")) +
          labs(x="Temperature treatment (°C)",y="Final germination proportions")

# View Figure 1: 'FGPfig'
FGPfig
# Save and export Figure 1: 'FGPfig' as a '.tiff' file
# This will automatically be placed in your working directory
ggsave(filename = "Fig 1 Final germination proportions.tiff", plot = FGPfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)


# Step 2: Estimating germination rates from the cumulative germination curves

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#????????????????????????????????????????????????????????????????????????
# Test the fit of different cumulative germination functions to your germination curves
fits <- drm(G/PG ~ Time, data = data, type = "binomial", fct = LL.2())
# create a table of all fits and call it 'selected'
selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                               LL2.5(), AR.2(), MM.2(), MM.3()), 
                    nested = FALSE, sorted = c("IC", "Res var", 
                                               "Lack of fit", "no"), linreg = F, icfct = AIC)
# The model listed in the first row is the model of best fit
# Where indicated in the following steps you must change the model to the one that
# best fits your data
#View, 'selected'
selected # In this case, it is model, 'W1.4'
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#????????????????????????????????????????????????????????????????????????

# Write the function to repeat what you just did above and call it FS
FS <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                          LL2.5(), AR.2(), MM.2(), MM.3()), 
                          nested = FALSE, sorted = c("IC", "Res var", 
                          "Lack of fit", "no"), linreg = F, icfct = AIC)
      row.names(selected)[1]} 

FSfit<-ddply(data, .(Grouping, Treatment), failwith(f = FS, quiet = T))

# View the table 'FSfit', which tells you the best fit across all 'Groupings', & 'Treatments'
# In this case the model, 'LL.2' best fits everything
FSfit

# Plot cumulative germination and check function fit visually
# If using any dataset other than the example dataset, you must take 
# action here and change 'fct=LL.2()' to the one indicated in table 'FSfit'
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

# Create a figure plotting all of your germination curves by treatment for every treatment
# by grouping and call the figure 'CGfig'
CGfig <- ggplot(data, aes(x = Time1, y = G/PG)) +
         geom_point(size=0.8, alpha = 0.5)  +
         geom_line(data = CGfit, 
                   aes( x = demo.fits.conc, y = Prediction)) +
         geom_ribbon(data = CGfit, 
                     aes(x = demo.fits.conc, y = Prediction, 
                     ymin = Lower, ymax = Upper), 
                     alpha = 0.2) + 
         facet_grid(Treatment ~ Grouping) +
         theme_bw() + theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(colour="black", fill="white")) +
         labs(x="Time",y="Final germination proportions")

# Look at Figure 2: CGfig
CGfig
# Save Figure 2: CGfig to a '.tiff' file in your working directory
ggsave(filename = "Fig 2 Cumulative germination curves.tiff", plot = CGfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)


# Set a function to fit the cumulative germination function 
# As above in 'CGfun' here you must change 'fct=LL.2()' to the appropriate model 
# that fits your data as indicated in the table 'FSfit'
GRfun <- function(x) {
         fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         doses <- ED(fits, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F)
                      }

# Apply the function to each treatment and add decils and rates and call it GR
GR <- ddply(data, .(Grouping, Treatment), failwith(f = GRfun, quiet = T))

# Add a new column,'Times' to the new dataset 'GR'
colnames(GR)[3] <- "Times"

# Calculate the inverse of time to obtain the germination rate for each treatment
GR$Rates <- 1/GR$Times
# and for each decil
GR$Decils <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
                     
# View Table 2: 'GR' 
View(GR)
# Write and export Table 2 'GR' to a new file
write.table(GR, "Table 2 Germination rates.txt", sep = "\t", col.names = NA)



# Step 3a: Fitting a segmented model to the full germination temperature range 

# Set a function 'CTfun'  to fit segmented regression and calculate cardinal 
# temperatures and thermal time
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
                    Tc = - asupra / bsupra,
                    thetasub = 1/bsub,
                    thetasupra = -1/bsupra)
                     }

# Apply the function 'CTfun' across 'Grouping' and 'Decils' to calculate cardinal temperatures
# for each decil and call it 'CT'
CT <- ddply(GR, .(Grouping, Decils), failwith(f = CTfun, quiet = T))

# View Table 3: 'CT'
View(CT)
# Write and export Table 3 to a file
write.table(CT, "Table 3 Cardinal temperatures.txt", sep = "\t", col.names = NA)


# Set a function to calculate segmented regressions
BLfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         seg <- segmented(lin.mod, seg.Z = ~ Treatment, psi = 25)
         data.frame(Treatment = unique(x)$Treatment, 
                    Rates = broken.line(seg)$fit)
                     }
# Fit the function to each 'Grouping' and 'Decil'
BLfit <- ddply(GR, .(Grouping, Decils), failwith(f = BLfun, quiet = T))

# Create your third figure, and call it 'BLfig'
BLfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point(size=1) + 
         facet_grid(Grouping ~ Decils) +
         geom_line(data = BLfit, colour = "black") +
         theme_bw() + theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(colour="black", fill="white")) +
         labs(x="Temperature treatment (°C)",y="Germination rate")

# View Figure 3: BLfig 
BLfig
# Save Figure 3 to your working directory
ggsave(filename = "Fig 3 Germination rate vs T (segmented model).tiff", plot = BLfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)


#Step 3b: Fitting a linear model to the sub- or supraoptimal germination temperature range

# Set a function to fit segmented regression and calculate base or ceiling T
LMfun <- function(x) {
         lin.mod <- lm(Rates ~ Treatment, data = x)
         a <- lin.mod$coefficients[1]
         b <- lin.mod$coefficients[2]
         data.frame(TbORTc = - a / b,
                    theta = abs(1/b))
                     }

# Calculate cardinal temperatures for each decil
LM <- ddply(GR, .(Grouping, Decils), failwith(f = LMfun, quiet = T))

# View Table 4: LM
View(LM)
# Create and export Table 4
write.table(LM, "Table 4 Base or ceiling temperatures.txt", sep = "\t", col.names = NA)

# Plot the GR versus time for each decil
LMfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
         geom_point(size=0.7)  +
         geom_smooth(method = "lm", size = 0.7, se = T, color="black") +
         facet_grid(Grouping ~ Decils) +
         theme_bw() + theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(colour="black", fill="white")) +
         labs(x="Temperature treatment (°C)",y="Germination rate")

# View Figure 4: LMfig
LMfig
# Write and export Table 4 to a file
ggsave(filename = "Fig 4 Germination rate vs T (linear model).tiff", plot = LMfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)



# Data exploration only* 
# Plot a smooth line figure for irregular data
SMTHfig <- ggplot(GR, aes(x = Treatment, y = Rates)) +
  geom_point(size=0.7)  +
  geom_smooth( size = 0.7, se = T, color="black") +
  facet_grid(Grouping ~ Decils) +
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     strip.background = element_rect(colour="black", fill="white")) +
  labs(x="Temperature treatment (°C)",y="Germination rate")
# View the Figure SMTHfig
SMTHfig
# Save and export SMTHfig
ggsave(filename = "Fig 4 Germination rate vs T (smooth model).tiff", plot = SMTHfig,
       path = NULL, scale = 1, width = 173, height = 173, 
       units = "mm", dpi = 300)


