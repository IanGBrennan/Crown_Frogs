library(stringr)
library(metablastr)
library(BUSpaRse)

source("~/Documents/GitHub/Crown_Frogs/dropEmptyTaxa.R")
source("/Users/ianbrennan/My Drive/R.Analyses/Convenient Scripts/Condensing_Alignments.R")
source("~/Documents/GitHub/Crown_Frogs/combineOrthologues.R")

setwd("~/Documents/GitHub/Crown_Frogs/Microhylidae/Alignments")
dropEmptyTaxa(file.extension = ".fasta")

#setwd("~/Desktop/Myobatrachidae/Alignments")
#setwd("~/Desktop/Microhylidae/Alignments")
setwd("~/Desktop/MM_Alignments/Combined")
extract.targets(taxon="longest",
                filetype=".fasta",
                replace.missing=T,
                path="~/Desktop/MM_Alignments/Combined")


setwd("/Users/ianbrennan/Documents/GitHub/Crown_Frogs/Microhylidae")
cfiles <- dir(getwd(), pattern=".fasta")

head(cfiles)
for (k in 1:length(cfiles)){
  new.name <- strsplit(cfiles[[k]], "__")[[1]][2]
  name.call <- paste("mv", paste0(getwd(),"/",cfiles[[k]]), paste0(getwd(),"/",new.name))
  system(name.call)
}

setwd("/Users/ianbrennan/Desktop/Microhylid_Combo/MACSE_Alignments")
afiles <- dir(getwd(), pattern=".fasta")

nfiles <- setdiff(afiles, cfiles)
for (j in 1:length(nfiles)){
  move.name <- paste("cp", paste0(getwd(),"/",nfiles[[j]]), "/Users/ianbrennan/Documents/GitHub/Crown_Frogs/Microhylidae")
  system(move.name)
}


# Compare the Microhylidae and Myobatrachidae alignments
rbh <- blast_best_reciprocal_hit(
  query = "~/Documents/GitHub/Crown_Frogs/Microhylidae/Microhylidae_All_Loci.fasta",
  subject = "~/Documents/GitHub/Crown_Frogs/Myobatrachidae/Myobatrachidae_All_Loci.fasta",
  search = "nucleotide_to_nucleotide",
  task = "dc-megablast",
  output.path = tempdir(),
  db.import = F
)
write.csv(rbh, file="~/Documents/GitHub/Crown_Frogs/RBH_Microhylidae_Myobatrachidae.csv", row.names=F, quote=F)
rbh.fasta <- rbh; rbh.fasta$query_id <- paste0(rbh.fasta$query_id,".fasta"); rbh.fasta$subject_id <- paste0(rbh.fasta$subject_id,".fasta")
# Combine the Microhylidae and Myobatrachidae alignments
setwd("~/Desktop/MM_Alignments")
combineOrthologues(filetype=".fasta", blast.table=rbh, path="~/Desktop/MM_Alignments")

# Which Microhylid loci weren't represented in the Myobatrachids?
mi.files <- dir("~/Desktop/MM_Alignments", pattern="Microhylidae")
excluded.mi <- setdiff(mi.files, rbh.fasta$query_id)
# Now copy them into the alignments directory
for (k in 1:length(excluded.mi)){
  move.call <- paste("cp", paste0(getwd(),"/",excluded.mi[[k]]), paste0(getwd(),"/Combined"))
  system(move.call)
}

# Which Myobatrachid loci weren't represented in the Microhylids?
my.files <- dir("~/Desktop/MM_Alignments", pattern="Myobatrachidae")
excluded.my <- setdiff(my.files, rbh.fasta$subject_id)
# Now copy them into the alignments directory
for (k in 1:length(excluded.my)){
  move.call <- paste("cp", paste0(getwd(),"/",excluded.my[[k]]), paste0(getwd(),"/Combined"))
  system(move.call)
}

# Compare the MM alignments to Xenopus
xrbh <- blast_best_reciprocal_hit(
  query = "~/Desktop/MM_Alignments/Combined/MM_All_Loci.fasta",
  subject = "~/Documents/GitHub/Crown_Frogs/Xenopus_tropicalis_cds.fasta",
  search = "nucleotide_to_nucleotide",
  task = "dc-megablast",
  output.path = tempdir(),
  db.import = F
)
write.csv(xrbh, file="~/Documents/GitHub/Crown_Frogs/RBH_MM_Xenopus.csv", row.names=F, quote=F)
xrbh.fasta <- xrbh; xrbh.fasta$query_id <- paste0(xrbh.fasta$query_id,".fasta"); xrbh.fasta$subject_id <- paste0(xrbh.fasta$subject_id,".fasta")
names(xrbh)[[1]] <- "subject_id"
names(xrbh)[[2]] <- "query_id"



# rename the Xenopus loci and extract just the ones that match our AHE targets
xt <- read.FASTA("~/Documents/GitHub/Crown_Frogs/Xenopus_tropicalis_cds.fasta")
locus.names <- names(xt)
names(xt) <- sapply(locus.names, function(x) str_split(x, " ")[[1]][1])
matches <- xt[which(names(xt) %in% xrbh$subject_id)]
write.FASTA(matches, file="~/Documents/GitHub/Crown_Frogs/Xenopus_tropicalis_AHE.fasta")

xn <- read.FASTA("~/Documents/GitHub/Crown_Frogs/Xenopus_tropicalis_AHE.fasta")
for (j in 1:length(xn)){
  locus.name <- names(xn)[[j]]
  align <- xn[j];
  names(align) <- "Xenopus_tropicalis_genome"
  write.FASTA(align, file=paste0("~/Desktop/MM_Alignments/Xenopus/",locus.name,".fasta"))
}

# Combine our Xenopus sequences with the MM_Alignments
setwd("~/Desktop/MM_Alignments/Combined")
combineOrthologues(filetype=".fasta", blast.table=xrbh, path="~/Desktop/MM_Alignments/Combined")

# Which MM loci weren't represented in the Xenopus?
mm.files <- dir("~/Desktop/MM_Alignments/Combined", pattern=".fasta")
excluded.mm <- setdiff(mm.files, xrbh.fasta$query_id)
# Now copy them into the alignments directory
for (k in 1:length(excluded.mm)){
  move.call <- paste("cp", paste0(getwd(),"/",excluded.mm[[k]]), paste0(getwd(),"/Combined_Xenopus"))
  system(move.call)
}

setwd("~/Desktop/MM_Alignments/Combined/Combined_Xenopus")
bads <- dir(getwd(), pattern="combined.")
for (y in 1:length(bads)){
  new.name <- paste0(str_split(bads[[y]], "__")[[1]][2],"__",str_split(bads[[y]], "__")[[1]][3])
  new.call <- paste("mv", paste0(getwd(),"/",bads[[y]]), paste0(getwd(),"/",new.name))
  system(new.call)
}
