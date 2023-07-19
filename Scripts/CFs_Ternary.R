library(Ternary)
library(patchwork)
library(RColorBrewer)
setwd("~/Crown_Frogs")

# read in the CF file and appropriate tree (ASTRAL)
support = read.delim("Trees/sCF/sCF_gCF.cf.stat", header=T, comment.char="#") # astral tree
astral.tree <- read.tree("Trees/Crown_Oz_ASTRAL.tre")


# extract support values for each clade
support.pelodryadidae <- dplyr::filter(support, ID >= 174 & ID <= 197)
support.myobatrachidae <- dplyr::filter(support, ID >= 141 & ID <= 164)
support.microhylidae <- dplyr::filter(support, ID >= 115 & ID <= 139)
support.asterophryinae <- dplyr::filter(support, ID >= 124 & ID <= 139)


support <- support.myobatrachidae


#### What we want to do us use a ternary plot to look at the gCF scores
# but gCFs don't add up to 100 because some loci don't provide definitive support
# so we will scale the remaining loci

#### Ternary Plot of rescaled Gene Concordance Factors
# choose the site-concordance-factor columns
gCF <- dplyr::select(support, "gCF", "gDF1", "gDF2"); gCF <- gCF[2:nrow(gCF),]
rs.gCF <- data.frame(t(sapply(1:nrow(gCF), function(x) (gCF[x,]/rowSums(gCF)[[x]])*100)))
colnames(rs.gCF) <- colnames(gCF)
rs.gCF <- mutate_all(rs.gCF, function(x) as.numeric(as.character(x)))
# set the background colors for point density
plotCols <- colorRampPalette(brewer.pal(6, "RdYlBu"), alpha=T); point.colors <- rev(plotCols(100))
# plot the ternary background
TernaryPlot(alab = 'gCF', blab = 'gDF1', clab = 'gDF2')
#ColourTernary(TernaryDensity(sCF, resolution = 10))
# plot the ternary point density colors
ColourTernary(TernaryDensity(rs.gCF, resolution = 10), spectrum=point.colors)

# plot the points, but start with ones where gCF > gDF1 & gDF2
TernaryPoints(dplyr::filter(rs.gCF, gCF > gDF1 & gCF > gDF2), col = "white", pch = 16, cex=1)
# plot points where gCF < gDF1 & 2
TernaryPoints(dplyr::filter(rs.gCF, gCF < gDF1 | gCF < gDF2), col = "black", pch = 16, cex=1)
# plot points where gCF~gDF1~gDF2
TernaryPoints(dplyr::filter(rs.gCF, gCF < 40 & gCF > 30 &
                              gDF1 < 40 & gDF1 > 30 &
                              gDF2 < 40 & gDF2 > 30), col = "grey", pch = 16, cex=1)



#### Ternary Plot of Site Concordance Factors
# choose the site-concordance-factor columns
sCF <- dplyr::select(support, "sCF", "sDF1", "sDF2"); sCF <- sCF[2:nrow(sCF),]
# set the background colors for point density
plotCols <- colorRampPalette(brewer.pal(6, "RdYlBu"), alpha=T); point.colors <- rev(plotCols(100))
# plot the ternary background
TernaryPlot(alab = 'sCF', blab = 'sDF1', clab = 'sDF2')
#ColourTernary(TernaryDensity(sCF, resolution = 10))
# plot the ternary point density colors
ColourTernary(TernaryDensity(sCF, resolution = 10), spectrum=point.colors)

# plot the points, but start with ones where gCF > gDF1 & gDF2
TernaryPoints(dplyr::filter(sCF, sCF > sDF1 & sCF > sDF2), col = "white", pch = 16, cex=1)
# plot points where gCF < gDF1 & 2
TernaryPoints(dplyr::filter(sCF, sCF < sDF1 | sCF < sDF2), col = "black", pch = 16, cex=1)
# plot points where gCF~gDF1~gDF2
TernaryPoints(dplyr::filter(sCF, sCF < 40 & sCF > 30 &
                            sDF1 < 40 & sDF1 > 30 &
                            sDF2 < 40 & sDF2 > 30), col = "grey", pch = 16, cex=1)


# investigate which the nodes are where an alternate topology is supported by CFs
filter(support, sCF < sDF1 | sCF < sDF2)
# investigate which nodes are equally supported in 3 topologies
filter(support, sCF < 40 & sCF > 30 &
         sDF1 < 40 & sDF1 > 30 &
         sDF2 < 40 & sDF2 > 30)

nrow(filter(support, sDF1 > sCF & sDF2 > sCF))
nrow(filter(support, sDF1 > sCF & sDF2 < sCF))



#### Investigate Gene Concordance Factors
# we won't Ternary plot because gCF+gDF1+gDF2 doesn't = 100

# identify nodes with low gCF values
filter(support, gCF < 40)
# investigate which the nodes are where an alternate topology is supported by CFs
filter(support, gCF < gDF1 | gCF < gDF2)




source("Scripts/color.by.CF.R")
# read in the mcmctree
mcmc.tree <- read.nexus("Trees/CrownFrogs_mcmctree.tre")
mcmc.tree <- bind.tip(mcmc.tree, tip.label="Heleophryne_purcelli", where = 11, edge.length = 2)

# for ASTRAL
CF.obj <- color.by.CF(astral.tree, cf.file="Trees/sCF/sCF_gCF.cf.stat", value.check = T, legend=T, CF="gCF")
# for MCMCtree
CF.obj <- color.by.CF(mcmc.tree, cf.file="Trees/sCF/sCF_gCF.cf.stat", value.check = T, legend=T)


SxG <- CF_comparison(CF.obj, type="sCFxgCF")
SxBL <- CF_comparison(CF.obj, type="sCFxBL")
GxBL <- CF_comparison(CF.obj, type="gCFxBL")

SxG + SxBL + GxBL


pelo.tree <- extract.clade(mcmc.tree, node=174)
pelo.data <- dplyr::filter(combo.edges, parent >= 174 & parent <= 197)
plot.phylo(pelo.tree, edge.color=pelo.data$CF_col, edge.width=3, cex=label.size,
           type="fan", open.angle=170); # can add 'edge.lty=combo.edges$gCF_lty' 










par(mar = rep(0.2, 4))
TernaryPlot(alab = 'a', blab = 'b', clab = 'c')

nPoints <- 4000L
coordinates <- cbind(abs(rnorm(nPoints, 2, 3)),
                     abs(rnorm(nPoints, 1, 1.5)),
                     abs(rnorm(nPoints, 1, 0.5)))




ggplot(gCF.obj, aes(x=gCF, y=sCF, label=child, fill=BranchLength)) +
  geom_point(shape=21, size=10) + theme_classic() + geom_smooth(color="black",linetype="dashed", level=0) + 
  scale_fill_fermenter(palette="RdYlBu", direction=1) + xlim(0,100) + ylim(0,100) + theme(legend.position = "bottom")






