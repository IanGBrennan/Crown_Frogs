library(ggplot2)
library(ggrepel)
source("Scripts/ggplotRegression.R")

rad <- read.csv("Scripts/Comparative_Radiations.csv", h=T)

ggplot() +
  geom_smooth(data=rad, aes(x=CrownAge, y=Species), method="glm", color="black") +
  geom_point(data=rad, aes(x=CrownAge, y=Species, color=Kind), size=3) +
  theme_bw() + geom_label_repel(data=rad, aes(x=CrownAge, y=Species, label = Group),
                                box.padding   = 0.35, 
                                point.padding = 0.5,
                                segment.color = 'grey50',
                                max.overlaps = 20)

rad[order(-rad$CrownAge),]

cage <- lm(data=rad, Species~CrownAge)

ggplotRegression(cage)
