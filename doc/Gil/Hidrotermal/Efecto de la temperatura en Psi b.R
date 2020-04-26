library(tidyverse)

data.frame(T = c(14, 16, 18, 20, 22, 24, 27, 28), 
           Psi = c(-1.64, -1.56, -1.45, -1.39, -1.24, -1.06, -0.63, -0.49)) %>%
ggplot(aes(x = T, y = Psi)) +
geom_smooth(method = "loess", alpha = 0.2) +
geom_point(size = 3) +
geom_vline(xintercept = 19.3, linetype = "dashed", color = "indianred") +
geom_text(aes(label = "To =  19.3ºC", x = 20, y = -1.65, color = "indianred"), show.legend = F) +
labs(x = "Temperatura (ºC)", y = expression(Psi[b]*(MPa))) -> fig1

ggsave(filename = "Efecto de la temperatura en Psi b.png", fig1, path = NULL, 
       scale = 1, width = 297, height = 180, units = "mm", dpi = 600)
