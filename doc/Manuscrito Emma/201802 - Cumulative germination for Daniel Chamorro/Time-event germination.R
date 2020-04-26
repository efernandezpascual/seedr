library(drc); library(tidyverse)

read.csv("Datos.csv") -> 
Datos

## AJUSTAR EL MODELO 

drm(Germinated ~ Start + End, data = Datos, type = "event", fct = LL.2()) ->
Modelo

## CALCULAR EL T50

ED(Modelo, 50)

## PREPARAR LOS DATOS PARA EL PLOT

expand.grid(conc = exp(seq(log(1.00e-04), log(max(Datos$Start)), length = 100))) ->
df1

data.frame(predict(Modelo, newdata = df1, interval = "confidence")) -> 
df2
 
df2$demo.fits.conc <- df1$conc 

Datos %>%
group_by(Replicate) %>%
mutate(Start =  ifelse(Start == 0, 1.00e-09, Start),
       Cumulative = cumsum(Germinated)/sum(Germinated)) %>%
filter(End != Inf) ->
df3

## DIBUJAR EL PLOT

ggplot(df2, aes(x = demo.fits.conc, y = Prediction)) +
geom_point(data = df3, aes(x = (Start + End)/2, y =  Cumulative, color = Replicate), size = 0.8, alpha = 0.5)  +
geom_line() +
geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) + 
theme_bw() + 
theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour = "black", fill = "white")) +
scale_y_continuous(labels = scales::percent) +
coord_cartesian(ylim = c(0, 1)) +
labs(x = "Time (days)",y = "Germination")


