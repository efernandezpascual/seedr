library(dplyr); library(purrr)

read.table("event_data.txt", header = T) %>%
filter(Grouping == "A" & Treatment == 23.75) %>%
group_by(Grouping, Treatment) %>%
do(possibly(GRfun(.)))

# Mirar como poder usar el do y evitar los errores