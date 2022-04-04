source("~/Documents/GitHub/Crown_Frogs/dropEmptyTaxa.R")
source("~/Google.Drive/R.Analyses/Convenient Scripts/Condensing_Alignments.R")

setwd("~/Desktop/Crown_Frogs/Myobatrachidae/Original_Alignments")
dropEmptyTaxa(file.extension = ".fasta")

setwd("~/Desktop/Crown_Frogs/Myobatrachidae/NT_Alignments")
extract.targets(taxon="Uperoleia_laevigata_MM1227",
                filetype=".fasta",
                replace.missing=T,
                path="~/Desktop/Crown_Frogs/Myobatrachidae/NT_Alignments")
