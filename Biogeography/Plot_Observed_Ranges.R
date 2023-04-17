library(phytools)
library(RColorBrewer)

setwd("CrownFrogs/Biogeography")

# read tree
btree <- read.nexus("CrownFrogs_BGB_Fossils_rotated.tre")
# rotate it so microhylids are on the outside
btree <- rotate(btree, node=108)

# get all the tip naems
all.tips <- btree$tip.label
# drop the outgroups
btree <- drop.tip(btree, tip=all.tips[93:102])
# get rid of the pelo ancestor
btree <- drop.tip(btree, tip="Pelodryadidae_Ancestor")

# read in the observed ranges
ranges <- read.csv("CrownFrogs_BGB_H1_ForPlotting.csv", header=T, row.names = 1)

# plot the tree
plotTree(btree,ftype="i",offset=5,fsize=0.5, 
         type="fan", part=0.5)
ranges <- ranges[btree$tip.label,]
tiplabels(pie=ranges,piecol=brewer.pal(8, "Spectral"),cex=0.3, pt.cex=1)
legend("bottom",colnames(ranges),pch=21,pt.bg=brewer.pal(8, "Spectral"),
       pt.cex=1, horiz=T)
