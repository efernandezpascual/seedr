# Load libraries
library(drc)
library(plyr)
library(segmented)
library(tidyr)
library(ggplot2)

# Load data
setwd("/Users/emmaladouceur/Desktop")
data <- read.table("dataGnsup_EM.txt", header = T)

View(data)
# Test the fit of different cumulative germination functions
fs <- function(x) {
      fits <- drm(G/PG ~ Time, data = x, type = "binomial", fct = LL.2())
      selected <- mselect(fits, list(LL.2(), LL.5(), W1.4(), W2.2(), LL2.2(), 
                          LL2.5(), AR.2(), MM.2(), MM.3()), 
                          nested = FALSE, sorted = c("IC", "Res var", 
                          "Lack of fit", "no"), linreg = F, icfct = AIC)
      row.names(selected)[1]}

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

View(times)


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


#rates$treatments
#change name of colum std. error
colnames(times)[which(names(times) == "Std. Error")] <- "stderror"
 
#spread  times data into wide format
times_wide<- spread(times, decils, rates)
View(times_wide)

#plot rates vs. treatment for t50
t50r <- ggplot(times_wide, aes(fill=factor(Treatment), y= t50, x=factor(Treatment))) 
dodge <- position_dodge(width=0.9) 
t50r + geom_bar(stat="identity", position=dodge) + scale_fill_manual( values=c("#0C5BB0FF", "#EE0011FF", "#15983DFF", "#EC579AFF" ,"#FA6B09FF" ,"#149BEDFF", "#A1C720FF", "#FEC10BFF","#16A08CFF", "#9A703EFF"))


#calculate cumulative germination 
#make with plyr
View(data)

#add cumulative germination percentage in a new column
cgp_data<-ddply(data, .(Treatment), transform, cgp = G/PG)


#Plot cumulative germination against time
germ<-ggplot(cgp_data, aes(x = Time, y = cgp)) + ylim(-0.03, 1.00) 
germ + geom_smooth(aes(color=factor(Treatment)), 
               se=FALSE, fullrange=TRUE) + scale_color_manual( values=c("#0C5BB0FF", "#EE0011FF", "#15983DFF", "#EC579AFF" ,"#FA6B09FF" ,"#149BEDFF", "#A1C720FF", "#FEC10BFF","#16A08CFF", "#9A703EFF"))


#regressions?????

ggplotRegression <- function (fit) {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + geom_point() + stat_smooth(col = "#74A089") + labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "Intercept =",signif(fit$coef[[1]],5 ), " Slope =",signif(fit$coef[[2]], 5), " P =",signif(summary(fit)$coef[2,4], 5)))}
ggplotRegression(lm(t50 ~  Treatment  , data = times_wide))


lin.mod <- lm(t50~Treatment, data=times_wide)
segmented.mod <- segmented(lin.mod, seg.Z = ~Treatment, psi=20)
summary(segmented.mod)

#segmented
ggplot(times_wide, aes(x = Treatment, y = t50)) +
  geom_point() +
  geom_line(data = segmented.mod, color = 'blue')


qplot(Treatment, t50, group = Treatment > 27, geom = c('point', 'smooth'), 
      method = 'lm', se = F, data = times_wide)

