library(phytools)

setwd("/Users/ianbrennan/Documents/GitHub/Crown_Frogs")

# read in original tree
ftree <- read.tree("Biogeography/CrownFrogs_BGB.tre")

# bind a tip for Baurubatrachus pricei or Calyptocephalella satan
ftree <- bind.tip(ftree, tip.label="Calyptocephalella_satan", 
                  where=which(ftree$tip.label=="Calyptocephalella_gayi"),
                  position=78, edge.length=8)

# bind a tip for the Antarctic Calypto, sister to C.gayi
ftree <- bind.tip(ftree, tip.label="Calyptocephalella_sp", 
                  where=which(ftree$tip.label=="Calyptocephalella_gayi"),
                  position=50, edge.length=10)

# bind a tip for Calyptocephalella canqueli, sister to C.gayi
ftree <- bind.tip(ftree, tip.label="Calyptocephalella_canqueli", 
                  where=which(ftree$tip.label=="Calyptocephalella_gayi"),
                  position=23, edge.length=2)

write.tree(ftree, file="Biogeography/CrownFrogs_BGB_Fossils.tre")

# find the MCRA of the Pelodryadidae
pelo <- getMRCA(ftree, tip=c("Litoria_citropa",
                              "Litoria_tyleri"))
# bind a tip for the hypothetical Antarctic Pelodryadid
ftree <- bind.tip(ftree, tip.label="Pelodryadidae_Ancestor",
                  where=pelo, position=6, edge.length=0.0000001)

plot(ftree, cex=0.3); axisPhylo()

write.tree(ftree, file="Biogeography/CrownFrogs_BGB_Fossils.tre")
