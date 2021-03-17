# Take the 5UTR, CDS and 3UTR sequences  
# Stick them together 
# Produce a new fasta file that has sequences in the format 5UTR-CDS-3UTR

# inputs 
# S_pombe 5UTR fasta, CDS fasta, 3UTR fasta. 
# these should have sequences for individual genes/exons, labelled with IDs
# Downloaded from PomBase

# Process 
# check that sequences are labelled with IDs
# for loop to add 5UTR to CDS
    # for sequence i in CDS file
    # store the ID of sequence i
    # find the sequence with the matching ID in 5UTRs
         # some may not have 5UTRs, if there isn't then i ++
    # merge the sequences into the format 5UTR-CDS
    # store product in position i of new object 5UTR-CDS
    # return 5UTR-CDS
    # i ++
# repeat above for adding 3UTR to 5UTR-CDS
    # for sequence i in 5UTR-CDS file
    # store the ID of sequence i
    # find the sequence with the matching ID in 3UTRs
          # some may not have 5UTRs, if there isn't then i ++
    # merge the sequences into the format 5UTR-CDS-3UTR
    # store product in position i of new object 5UTR-CDS-3UTR
    # return 5UTR-CDS-3UTR
    # i ++

# output a fasta file with the UTRS and CDS for genes.  
library(Biostrings)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(parallel)
library(rhdf5)

three_UTR <- readDNAStringSet('3UTR.fa.gz')
CDS <- readDNAStringSet('cds.fa.gz')
five_UTR <- readDNAStringSet('5UTR.fa.gz')


CreateFastaFullUTRs <- function(three_UTR, CDS, five_UTR)
{
  cc <- 1
  fiveutr_CDS_threeutr <- list()
  # loop to add, if present, 5UTRs upstream of CDS sequences and 3UTRs downstream.
  # add names column 
  for (i in unique(names(CDS))){
    fiveutr_CDS_threeutr[[cc]] <- DNAStringSet(paste0(unlist(
        as.character(five_UTR[names(five_UTR)==i])),
        unlist(as.character(CDS[names(CDS)==i])), 
        unlist(as.character(three_UTR[names(three_UTR)==i]))
        ))
    cc <- cc + 1 
  }
  fiveutr_CDS_threeutr <- unlist(DNAStringSetList(fiveutr_CDS_threeutr))
  names(fiveutr_CDS_threeutr) <- unique(names(CDS))
  
  # check that the length of a sequence fiveUTR_CDS_seqlist is equal
  # to the combined lengths of seperate CDS and five_utr length for each gene 
  for (j in unique(names(fiveutr_CDS_threeutr))){
      if (width(fiveutr_CDS_threeutr[names(fiveutr_CDS_threeutr)==j])!= 
                                      (sum(width(five_UTR[names(five_UTR)==j])) +
                                          sum(width(CDS[names(CDS)==j])) + 
                                          sum(width(three_UTR[names(three_UTR)==j])))){
          print('Length of full sequence is not the sum of the length of the 3UTR, CDS and 5UTR')
          break
      } 
  }
  return(fiveutr_CDS_threeutr)
  print('all lengths okay')
}

# export fasta files with 3UTR, CDS and 5UTR sequences.  

writeXStringSet(fiveutr_CDS_threeutr,filepath = file.path('.','S_pombe_full_UTR.fasta'),format = "fasta")
