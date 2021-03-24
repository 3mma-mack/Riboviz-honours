# Take the 5UTR, CDS and 3UTR sequences  
# Stick them together 
# Produce a new fasta file that has sequences in the format 5UTR-CDS-3UTR

# inputs 
# S_pombe 5UTR fasta, CDS fasta, 3UTR fasta. 
# these should have sequences for individual genes/exons, labelled with IDs
# Downloaded from PomBase: https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/

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

# run the function.

CreateFastaFullUTRs(three_UTR, CDS, five_UTR)

# update names, so match with names used in gff file later.

names(fiveutr_CDS_threeutr) <- paste(names(fiveutr_CDS_threeutr), '.1', sep = '')
# export fasta files with 3UTR, CDS and 5UTR sequences.  

writeXStringSet(fiveutr_CDS_threeutr,filepath = file.path('.','S_pombe_full_UTR.fasta'),format = "fasta")

# Make GFF ranges match UTR
# gff, 5UTR, CDS, 3UTR.
# make names consistent 
# for each gene, if type == UTR5, update the range to be 1 to length of 5UTR for that gene
# if type == CDS, update range to be length of 5UTR +1 to length of 5UTR + CDS 
# if type == UTR3, update range to be  ength 5UTR+CDS+1 to length of 5UTR+CDS+3UTR
# return updated gff 

genome <- readDNAStringSet('S_pombe_full_UTR.fasta',format = "fasta")
annot <- readGFFAsGRanges('Alex_files/Schizosaccharomyces_pombe_CDS_w_250utrs.gff3')

# make the names compatible to allow manipulation



names(five_UTR) <- paste(names(five_UTR), '.1', sep = '')
names(CDS) <- paste(names(CDS), '.1', sep = '')
names(three_UTR) <- paste(names(three_UTR), '.1', sep = '')

# seperate out the UTRs and CDS to allow processing 

gff_fiveutr <- annot[annot$type == 'UTR5']
gff_cds <- annot[annot$type == 'CDS']
gff_threeutr <- annot[annot$type == 'UTR3']

# Edit the ranges of the 5UTR gff to be 1-width of listed UTR for each gene, if present

  for(i in gff_fiveutr$Name){
    if(i %in% names(five_UTR)){
      ranges(gff_fiveutr[gff_fiveutr$Name == i]) <- IRanges(start = 1, 
                                                          width = width(five_UTR[names(five_UTR)==i]))
      }else{
         start(ranges(gff_fiveutr[gff_fiveutr$Name == i])) <- 0
         width(ranges(gff_fiveutr[gff_fiveutr$Name == i])) <- 1
         }
    }

# repeat for CDS, checking if a 5' UTR is present 

  for(i in gff_cds$Name){
    if(i %in% names(CDS)){
      if(i %in% names(five_UTR)){
        ranges(gff_cds[gff_cds$Name == i]) <- IRanges(start = end(ranges(gff_fiveutr[gff_fiveutr$Name == i])) +1,
                                                     width = width(CDS[names(CDS)==i]))  
         }else{
           ranges(gff_cds[gff_cds$Name == i]) <- IRanges(start = 1, width = width(CDS[names(CDS)==i]))  
            }
       }
    }

  for(i in gff_threeutr$Name){
    if(i %in% names(three_UTR)){
      ranges(gff_threeutr[gff_threeutr$Name == i]) <- IRanges(start = end(ranges(gff_cds[gff_cds$Name == i])) +1,
                                                            width = width(three_UTR[names(three_UTR)==i]))
       }else{
          ranges(gff_threeutr[gff_threeutr$Name == i]) <- IRanges(start = end(ranges(gff_cds[gff_cds$Name == i])), 
                                                            width = 1)
          }
    }

# combine GRangfes objects into one to be in the order 5utr - cds - 3utr

c <- 1
total_gff <- GRanges()

for(i in 1:length(gff_cds)){
  total_gff[c]<- gff_fiveutr[i]
  total_gff[(c+1)] <- gff_cds[i]
  total_gff[(c+2)] <- gff_threeutr[i]
  c <- c+3
}

#extra genes are present in annot compared to CDS, all of which are alternative transcripts for other genes. there are 8 in total, so remove
annot_names <- unique(annot$Name)
genome_names <- names(genome)
extra_genes <- annot_names[!(annot_names %in% genome_names)]
test <- Pasted_gff[-c(which(Pasted_gff$Name %in% c("SPAC17G6.02c.2", "SPBC2D10.10c.3", "SPBC2D10.10c.2", "SPCC162.04c.2",  "SPCC1620.02.2", "SPCC1906.03.2",  "SPCC548.03c.2" , "SPNCRNA.103.2" )))]

# export gff 
export.gff3(test, con=file.path('.','S_pombe_full_UTR_removed_genes_pasted.gff3'))

# as a quicker way to get final gff, paste together to be all the 5utr, then all CDS sequences, then 3utr
Pasted_gff <- GRanges(c(gff_fiveutr, gff_cds, gff_threeutr))
export.gff3(Pasted_gff, con=file.path('.','S_pombe_pasted_full_UTR.gff3'))

