library(tidyverse)

centaury %>%
  ggplot(aes(x =  Time, y = Germinated/Germinable, group = Temperature)) +
  facet_wrap(~ Population) +
  geom_line()

centaury %>% filter(Time <= 28) %>%
  filter(Population %in% c("CM", "MA")) %>%
  ggplot(aes(x =  Time, y = Germinated/Germinable, group = Temperature)) +
  facet_wrap(~ Population) +
  geom_line()


centaury %>% filter(Time <= 28) %>%
  filter(Population %in% c("CM", "MA")) %>%
  group_by(Population, Temperature) %>%
  mutate(Germinate = c(diff(c(0, Germinated)))) %>% write.csv("cent.csv")

diff(centaury$Germinated)
