library(metablastr)

setwd("Crown_Frogs/Alignments/Crown_Oz_Alignments/")
source("../../Scripts/Condensing_Alignments.R")

# make a fasta including all the AHE loci
extract.targets(taxon="longest", filetype=".fasta", replace.missing=T, path=".")

# 

# Compare the Microhylidae and Myobatrachidae alignments
rbh <- blast_best_reciprocal_hit(
  query = "Crown_Frogs/Alignments/Crown_Oz_Alignments/longest_All_Loci.fasta",
  subject = "Crown_Frogs/Alignments/Xenopus_tropicalis_exons.txt",
  search = "nucleotide_to_nucleotide",
  task = "dc-megablast",
  output.path = tempdir(),
  db.import = F
)

exon.hits <- rbh$subject_id
exon.hits <- sapply(exon.hits, function(x) strsplit(x, "\\|")[[1]][5])
names(exon.hits) <- NULL
rbh$subject_id_name <- exon.hits

write.csv(rbh, file="../RBH_AHE_XenopusExons.csv")


# Instructions for downloading the exons:
# 1. Go to BioMART website. (http://www.ensembl.org/biomart/martview/b3c3cbdd30499bc130975ee760be3b7c)
# 2. Choose database: Ensembl Genes 70
# 3. Choose dataset: Homo sapiens
# 4. Click 'Attributes' then select the 'Sequences' option
# 5. Expand the sequences pane and select the 'Exon sequences' option
# 6. Expand the 'Header information' pane. Select the info you want to associate with each sequence record (e.g. Ensembl Gene ID, Ensembl Transcript ID, Ensembl Exon ID).
# 7. Click the 'Results' button.
# 8. Select the 'Compressed File' option and hit the 'Go' button.