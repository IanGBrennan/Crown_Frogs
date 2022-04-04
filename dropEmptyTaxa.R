require(ips)

dropEmptyTaxa <- function(file.extension = ".fasta"){
  # read in the names for all the alignment files
  files <- dir(getwd(), pattern=file.extension)
  
  for (p in 1:length(files)){
    
    # get the locus name
    locus.name <- stringr::str_split(files[p], pattern = file.extension)[[1]][1]
    
    # read in the alignment file
    aligned <- read.dna(files[p], format="fasta", as.matrix=TRUE)
    
    # remove empty rows and columns
    #nonempty <- deleteEmptyTaxa(aligned, quiet=T)
    nonempty <- ips::deleteEmptyCells(aligned, margin=1, nset=c("-","n","?","N"))
    
    # make sure there's still sequences:
    if (length(nonempty) > 0) {
      
      # write the file (appending each new sequence)
      write.FASTA(nonempty, paste0("Reduced_", files[p]))
      
      # make a file which tells us which samples are missing from which alignments
      missing.msg <- paste(files[p], "lacks data for taxa:", setdiff(rownames(aligned), rownames(nonempty)))
      write.table(missing.msg, file=paste0("Missing_Samples.txt"), append=T, row.names=F, col.names=F, quote=F)
      
      # I started desiging a partition file, but easier to do with 'concat' in AMAS
      #part.info <- paste("charset", locus.name, "=", )
      
    } else {
      print(paste("alignment", files[p], "had no sequence data for any taxa"))
      missing.locus <- paste(files[p], "is now empty")
      write.table(missing.locus, file=paste0("Missing_Loci.txt"), append=T, row.names=F, col.names=F, quote=F)
    }
    
    print(paste(files[p],";", (p+1)))
  }
}



