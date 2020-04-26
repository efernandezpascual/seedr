# Open data in cumulative germination format
df1 <- read.table("S1 Example Data.txt", header = T)
df1 <- subset(df1, Grouping == "B" & Treatment == 18.75)[, 4:6]

# Transform to time-event format
df2 <- data.frame(Start = c(df1$Time), End = c(subset(df1, Time > 0)$Time, Inf),
                  Germinated = c(diff(df1$G), 
                                 max(df1$PG) - tail(df1$G, 1)),
                  Total = max(df1$PG))

m1 <- drm(G/PG ~ Time, data = df1, type = "binomial", fct = LL.2())
summary(m1)
ED(m1, 50)

m2 <- drm(Germinated ~ Start + End, data = df2, type = "event", fct = LL.2())
summary(m2)
ED(m2, 50)