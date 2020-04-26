
#install these packages
library(effects)
library(MASS)
library(plyr)

#set your own working directory
setwd("/Users/emmaladouceur/Desktop/MUSE/A Research/Ob 5/scripts")

#set your own dataset
data<-read.table("dataGnsup_EM.txt", header = T)

#p -> proportion of germination, dose-> time
times <- dlply(data, .(Treatment), function(data) {
          GLM <- glm(cbind(G, PG - G) ~ Time, data, family=binomial,control = list(maxit = 50))
          dose.p(GLM, p = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
          })
times


##Germination rate (1/t50)- slope, graph

rates <- ddply(data, .(Treatment), function(data) {
          GLM <- glm(cbind(G, PG - G) ~ Time, data, family=binomial,control = list(maxit = 50))
          t <- dose.p(GLM, p = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
          1/t
          })
rates
#change to your working directory
write.csv(rates, "/Users/emmaladouceur/Desktop/rates.csv" )




#########********MORE COMPICATED STUFF******############
## change time, model with final germination (end of experiment on day "23" only)
m1 <- glm(cbind(G, PG - G) ~ Treatment, subset(data, Time == "30"), family=binomial,control = list(maxit = 50))
summary(m1)

#excel draw (effect size)
#cut test another fresh and dead seeds
proportions <- data.frame(Effect(c("Treatment"), m1)$x, 
               Est. = plogis(Effect(c("Treatment"), m1)$fit), 
               lower = plogis(Effect(c("Treatment"), m1)$lower),
               upper = plogis(Effect(c("Treatment"), m1)$upper))
proportions


plot(proportions)


#install this if you do not have it
library(ggplot2)

#you can tell R what order you want to see the treatments in
#(lowest-highest) fill this in yourself.
proportions$Treatment <- factor(proportions$Treatment, levels=c("2C","7-2C","10-5C","15-5C","17-8C","20C","25-15C","30C","40C"))

#make your effect plots based on the 'propertions' we just calculated


###***CHECK OUT THIS PLOT
#the whiskers are not confidence intervas, but the lower and upper estimate of the
#effect prediction
limits <- aes(ymax = proportions$upper, ymin=proportions$lower)
pa <- ggplot(proportions, aes(fill=Treatment, y=Est., x=Treatment)) + ylim(0.00, 1.00) 
dodge <- position_dodge(width=0.9)
pa + geom_bar(stat="identity",position=dodge) + geom_errorbar(limits,  position=dodge, width=0.25, lwd=2) +scale_fill_manual(values=c("#0C5BB0FF", "#EE0011FF", "#15983DFF", "#EC579AFF" ,"#FA6B09FF" ,"#149BEDFF", "#A1C720FF", "#FEC10BFF","#16A08CFF", "#9A703EFF"))




#LAST STEP
#LINEAR REGRESSION
#change the columns headings in the 'rates' file you created to just be numbers

datar<-read.csv("rates.csv", header = T)
View(datar)


#tell r which order to put them in
#order them from lowest-highest temp
datar$Treatment <- factor(datar$Treatment, levels=c("2","7","10","15","17-8C","20C","25-15C","30C","40C"))


#CURVED REGRESSION
ggplotRegression <- function (fit) {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + geom_point() + stat_smooth(col = "#74A089") + labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), "Intercept =",signif(fit$coef[[1]],5 ), " Slope =",signif(fit$coef[[2]], 5), " P =",signif(summary(fit)$coef[2,4], 5)))}
ggplotRegression(lm(t50 ~  Treatment  , data = datar))

#use segmented package to create segments
library(segmented)
lin.mod <- lm(t50~Treatment, data=datar)
segmented.mod <- segmented(lin.mod, seg.Z = ~Treatment, psi=20)
summary(segmented.mod)

#segmented
ggplot(datar, aes(x = Treatment, y = t50)) +
  geom_point() +
  geom_line(data = segmented.mod, color = 'blue')


#create my own break
qplot(Treatment, t50, group = Treatment > 25, geom = c('point', 'smooth'), 
      method = 'lm', se = F, data = datar)



