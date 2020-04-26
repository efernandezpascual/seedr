
# An easy and automated calculation of the germination cardinal temperatures and thermal time using R
# Authors; Emma Ladouceur, Hugh Pritchard, and Eduardo Fernandez-Pascual

# Beginners please refer to Supplementary File 3 for support in getting started

# Install packages & dependencies
install.packages("binom", dependencies=TRUE)
install.packages("dplyr", dependencies=TRUE)
install.packages("drc", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("plyr", dependencies=TRUE)
install.packages("segmented", dependencies=TRUE)

# Load libraries
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
data <- read.table("Supplementary Information  1_29_12_EL.txt", header = T)

# Explore the data
  # View the data
View(data) #or
data
  # How many columns and rows does the dataset have?
ncol(data) #columns
nrow(data) #rows
  # What are the column names?
colnames(data)
  # What are the 'levels' of the column, 'Grouping' in the dataset 'data' ?
levels(data$Grouping)


# Ready to begin!

# Step 1: Checking whether the data represents the full germination temperature range

# First filter data by the final scoring date
FSD <- data %>% 
       group_by(Treatment, Dish) %>%
       filter(Time == max(Time))

# Write a function to estimate mean germination proportions and the associated confidence interval
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
# this will automatically be written to the working directory you already sent 
write.table(FGP, "Table 1 Final germination proportions.txt", sep = "\t", col.names = NA)

# Prepare to plot final germination proportions

# Format all future plots
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

# Test the fit of different cumulative germination functions to your germination curves

FS <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                          LL2.5(), AR.2(), MM.2(), MM.3()), 
                          nested = FALSE, sorted = c("IC", "Res var", 
                          "Lack of fit", "no"), linreg = F, icfct = AIC)
      row.names(selected)[1]} 


FSfit<-ddply(data, .(Grouping, Treatment), failwith(f = FS, quiet = T))

# Plot cumulative germination and check function fit visually
# If using any dataset other than the example dataset, take 
#action here and change 'fct=LL.2()' to.....
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
# that fits your data
GRfun <- function(x) {
         fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
         doses <- ED(fits, c(10, 20, 30, 40, 50, 60, 70, 80, 90), display = F)
                      }

# Apply the function to each treatment and add decils and rates and call it GR
GR <- ddply(data, .(Grouping, Treatment), failwith(f = GRfun, quiet = T))

# Add a column, 'Times' to the new dataset 'GR'
colnames(GR)[3] <- "Times"

# Calculate the inverse of time to obtain the germination rate for each treatment
GR$Rates <- 1/GR$Times
#  and for each decil
GR$Decils <- c("t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90")
                     
# View Table 2: 'GR' 
View(GR)
# Write and export  Table 2 'GR' to a new file
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
# Smooth line figure for irregular data
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


