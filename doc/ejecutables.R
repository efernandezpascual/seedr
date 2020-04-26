devtools::document()
#usethis::use_vignette("seedr_Hydro_and_Thermal_Time_Germination_Models")
devtools::build()
# usethis::use_data(prueba, prueba, overwrite = T)

# Physiodata data object

physiodata(centaury) -> cent1
physiodata(centaury, groups = "population") -> cent2
physiodata(grasses, x = "psi") -> grass1
physiodata(grasses, x = "psi", groups = "species") -> grass2

summary(cent1)
summary(cent2)
summary(grass1)
summary(grass2)

cent1
cent2
grass1
grass2

plot(cent1)
plot(cent2)
plot(grass1)
plot(grass2)

barplot(cent1)
barplot(cent2)
barplot(grass1)
barplot(grass2)

# Testing Bradford's model

physiotime(d = grasses, t = "times", g = "germinated", pg = "germinable", x = "psi",
           groups = c("species", "temperature"), method = "bradford") -> m1 # Dos grupos

physiotime(d = subset(grasses, species == "Anisantha rubens"), t = "times", g = "germinated", pg = "germinable", x = "psi",
           method = "bradford") -> m3 # Un grupo


m1 # print
m3


summary(m1) # proporciones de germinación por tratamiento
summary(m3)

plot(m1) # plot de los valores transformados según el modelo
plot(m3)



# Huidobro

physiotime(centaury, x = "temperature", method = "huidobro", groups = "population") -> mTemp
mTemp
summary(mTemp)
plot(mTemp)

physiotime(subset(centaury, population == "Malconcecho"), x = "temperature", method = "huidobro") -> mTemp1
mTemp1
summary(mTemp1)
plot(mTemp1)

physiotime(grasses, x = "psi", method = "huidobro", groups = c("species", "temperature")) -> mTemp2
mTemp2
summary(mTemp2)
plot(mTemp2)

# Para Gil

# Huidobro da un error en algunos casos que sólo hay suboptimo
physiodata(subset(centaury, temperature < 15), x = "temperature") -> cent1
huidobro(cent1$proportions) # Max R2 falla
huidobro(cent1$proportions, tops = "Max value") # Max value no

# El ajuste de Huidobro parece muy malo
# Tanto con datos de water potential (grasses) como de temperatura(centaury)
physiodata(centaury, x = "temperature") -> cent
physiodata(grasses, x = "psi") -> grass
huidobro(grass$proportions) -> m1
huidobro(cent$proportions) -> m2
m1
m2
plot(m1)
plot(m2)
# Bradford parece mejor para ambos casos
bradford(grass$proportions) -> m1
bradford(cent1$proportions) -> m2
m1
m2
plot(m1)
plot(m2)
