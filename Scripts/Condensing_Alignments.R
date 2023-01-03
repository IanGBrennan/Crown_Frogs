# in your terminal, delete the last line of each file, if there's an extra space/gap there with:
# $ sed -i '$d' PATH_TO_FILE
# or $ sed -i '$d' *.fasta # to do all fasta files

#### We'll make a loop to change the file format from Phylip to Fasta first
#path = "/Users/ianbrennan/Desktop/Microhylid_Combo"
#out.file<-""
#file.names <- dir(path, pattern =".fasta")
#dir.create(paste0(path,"/fasta")) # make a directory for the new files

extract.targets <- function(taxon, filetype, replace.missing=T,
                            path=NULL){
  file.names <- dir(path, pattern=filetype)
  taxon.og <- taxon
  for (p in 1:length(file.names)) {
    # first remove the first line with the number of taxa and characters
    file <- paste(path, "/", file.names[p], sep="")
    short.file <- file.names[p]
    shortie <- strsplit(short.file, filetype)[[1]][1]
    #removed <- paste("sed -i '1d'", file)
    #system(removed)
    
    # read in the alignment 
    if(filetype==".fasta"){full_alignment <- read.dna(file, format="fasta", as.matrix=TRUE)}
    else if(filetype==".phylip"){#taxon <- paste0(taxon.og,"\t"); 
                                 full_alignment <- read.dna(file, format="sequential")}
    
    if(taxon == "longest"){
      # get the number of missing bases per sample
      missing.sum <- apply(full_alignment, MARGIN = 1, FUN = function(x) length(which(as.numeric(x) %in% c(4, 240, 2))))
      new.taxon <- names(missing.sum[which(missing.sum == min(missing.sum))])
      if(length(new.taxon > 1)){new.taxon <- new.taxon[[1]]}
      
      # pull out just the new target taxon
      target.alignment <- full_alignment[new.taxon,]
      
      # rename the alignment with the locus name
      row.names(target.alignment) <- shortie
      
      # write the file (appending each new sequence)
      if (p==1){
        if(filetype==".fasta") {write.FASTA(target.alignment, paste0(path, "/", taxon, "_All_Loci.fasta"))}
        else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(path, "/", taxon.og, "_All_Loci.fasta"))}
      }
      else {
        if(filetype==".fasta") {write.FASTA(target.alignment, paste0(path, "/", taxon, "_All_Loci.fasta"), append=T)}
        else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(path, "/", taxon.og, "_All_Loci.fasta"), append=T)}
      }
    }
    
    else if (!taxon == "longest"){
      if(!taxon %in% row.names(full_alignment)){
        if(replace.missing==F){
          # tell us that the taxon is missing from the current alignment
          missing <- paste(taxon, "is not in alignment", shortie)
          print(missing)
        }
        else if(replace.missing==T){
          new.taxon <- rownames(full_alignment)[[1]]
          
          # tell us that the taxon is missing from the current alignment, and what we'll use instead
          missing <- paste(taxon, "is not in alignment", shortie, "replacing with sequence from", new.taxon)
          print(missing)
          
          # pull out just the new target taxon
          target.alignment <- full_alignment[new.taxon,]
          
          # rename the alignment with the locus name
          row.names(target.alignment) <- shortie
          
          # write the file (appending each new sequence)
          if(filetype==".fasta") {write.FASTA(target.alignment, paste0(path, "/", taxon, "_All_Loci.fasta"), append=T)}
          else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(path, "/", taxon.og, "_All_Loci.fasta"), append=T)}
        }
        # write the info to file
        write.table(missing, file=paste0(path,"/",taxon,"_MISSING_LOCI.txt"), append=T, row.names=F, col.names=F)
      }
      else if(taxon %in% row.names(full_alignment)){
        # pull out just the target taxon
        target.alignment <- full_alignment[taxon,]
        
        # rename the alignment with the locus name
        row.names(target.alignment) <- shortie
        
        # write the file (appending each new sequence)
        if(filetype==".fasta") {write.FASTA(target.alignment, paste0(path, "/", taxon, "_All_Loci.fasta"), append=T)}
        else if(filetype==".phylip") {ips::write.fas(target.alignment, paste0(path, "/", taxon.og, "_All_Loci.fasta"), append=T)}
      }
    }
  }
}

#extract.targets(taxon="Austrochaperina_fryi_1", filetype=".fasta")


 