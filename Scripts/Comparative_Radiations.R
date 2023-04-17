library(ggplot2)
library(ggrepel)
library(dplyr)
source("Scripts/ggplotRegression.R")

# read in the radiation information
rad <- read.csv("Scripts/Comparative_Radiations.csv", h=T)

# plot the relationship between crown age and species richness
ggplot() +
  #geom_smooth(data=rad, aes(x=CrownAge, y=Species), method="glm", color="black") +
  geom_point(data=rad, aes(x=CrownAge, y=Species, color=Kind), size=3) +
  theme_bw() + geom_label_repel(data=rad, aes(x=CrownAge, y=Species, label = Group),
                                box.padding   = 0.35, 
                                point.padding = 0.5,
                                segment.color = 'grey50',
                                max.overlaps = 20)

# fit a linear model to the data (excluding Limno and Myo, which are included in Myobatrachoidea)
cage <- lm(data=filter(rad, !Group %in% c("Limnodynastidae", "Myobatrachidae")), 
           Species~CrownAge)
# plot the linear model
ggplotRegression(cage)

# fit a linear model to just the frog data
fage <- lm(data=filter(rad, Group %in% c("Cophixalus", "Austrochaperina", "Ranidae",
                                         "Pelodryadidae", "Myobatrachoidea")),
           Species~CrownAge)

# plot the linear model
ggplotRegression(fage)

# or try it with pipes
rad %>%
  select(Species, CrownAge, Group) %>%
  filter(Group %in% c("Cophixalus", "Austrochaperina", "Ranidae",
               "Pelodryadidae", "Myobatrachoidea")) %>%
  lm(data=., Species~CrownAge) %>%
  ggplotRegression()
