require(ape)

#path = "/Users/Ian/Desktop/Microhylid_Combo" # make sure to designate you path appropriately!
#file.names <- dir(path, pattern =".fasta")
#setwd(path)


############## THIS IS THE NEW VERSION (FAST) ###################
#### It works by adding samples to an existing alignment, not building a new one together
#### Start by combining the loci that are orthologous, into a single Phylip file

combineOrthologues <- function(filetype = c(".fasta"), blast.table, path=NULL){
  file.names <- dir(path, pattern =".fasta")
  
  for (z in 1:nrow(blast.table)) {
    
    # select the locus names
    locus.a <- blast.table[z,"query_id"]
    locus.b <- blast.table[z,"subject_id"]
    
    # identify the files
    file.a <- paste(locus.a, filetype, sep="")
    file.b <- paste(locus.b, filetype, sep="")
    
    # determine if any of the query sequences have to be reverse complemented
    if(blast.table$q_start[[z]] > blast.table$q_end[[z]]){
      reverse.it <- paste("locus", locus.a, "has been reversed"); print(preverse.it)
      locus.in <- read.FASTA(file.a)
      locus.out <- rev(ape::complement(locus.in))
      file.a <- paste0(locus.a, "_RC", filetype)
      write.FASTA(locus.out, paste0("RC_",file.a))
      
      # document if any of the sequences have been flipped around
      write.table(reverse.it, file=paste0(path,"/","Combined_Alignments_RCs.txt"), append=T, row.names=F, col.names=F)
    }
    # determine if any of the subject sequences have to be reverse complemented
    if(blast.table$s_start[[z]] > blast.table$s_end[[z]]){
      reverse.it <- paste("locus", locus.b, "has been reversed"); print(reverse.it)
      locus.in <- read.FASTA(file.b)
      locus.out <- rev(ape::complement(locus.in))
      file.b <- paste0(locus.b, "_RC", filetype)
      write.FASTA(locus.out, file.b)
      
      # document if any of the sequences have been flipped around
      write.table(reverse.it, file=paste0(path,"/","Combined_Alignments_RCs.txt"), append=T, row.names=F, col.names=F)
    }
    
    new.name <- paste("combined.",locus.a,"__",locus.b,".fasta", sep="")
    
    call <- paste(paste0(getwd(), "/muscle3.8.31_i86darwin64"),
                  "-profile -in1 ", file.a, " -in2 ", file.b, " -out ", new.name)
    system(call)
  }
}

#combineOrthologues(filetype=".fasta", blast.table = rbh)
