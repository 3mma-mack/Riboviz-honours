# rename files to have same name, return files. general function to be applied to all

# This is dependent on the fasta and gff having a consistent style of labeling 
# ie 'chr II' will be relabelled to 'II' if the annot name is 'II' but 
# 'chr 2' will not be relabelled to 'II' 

# input: GRanges annot, DNAStringSet genome 

# find unique names from annot - in annot@seqnames@values - these are chromosome names
# Save unique names as new.names
# check that unique names match with part of genome names ie 'I' with 'chr I S.pombe'
   # loop; for each chromosome name (names(genome)), check that one of the new names 
   # matches with a seqment of the name  
        # rename genome name to with new names, where matching
   # if no matches return error and 'naming style are inconsistent,
   # may need to rename them prior to run'
   # return(renamed genome DNAStringSet)

# output: DNAStringSet genome


library(Biostrings)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(parallel)
library(rhdf5)


input <- "Schizosaccharomyces_pombe_all_chromosomes.fa.gz"
gff <- "Schizosaccharomyces_pombe_all_chromosomes.gff3.gz"

annot <- readGFFAsGRanges(gff)
genome <- readDNAStringSet(input,format = "fasta")

# rename for testing with different styles to see if the function works correctly  
names(genome) <- c("chr_II_telomeric_gap dd","CHR I Spombe"," CHR II Spombe"," CHR III Spombe","zz mating_type_region","s mitochondrial")
names(genome) <- c("chr_II_telomeric_gap dd","CHR 1 Spombe"," CHR 2 Spombe"," CHR III Spombe","zz mating_type_region","s mitochondrial")

# find the chromosome labels used in the gff file and extract those strings from the names of the genome
# The list new.names is made based on the order the names appear in the names(genome), not in the gff, so can just be returned following processing 
# If the naming style is inconsistent, ie '1' and 'I', this is shown as 'character(0)' in new.names
# in 'new.names', 'II' is listed as 'I''I'. collapsing to resolve
# resolving in a loop, producing a list that is only as long as the previous list of names
# using an empty list resulted in a length of 7 with [7] being filled with ''
# add the new names

renameGenome <- function(genome, annot)
{
   new.names <- str_extract_all(names(genome), paste(unique(as.character(annot@seqnames@values)),collapse="|"))
     if(any(new.names == 'character(0)')){
      print('Chromosome naming style may not be consistent, manually rename using names(genome) <- x')
      return(0)
     }
   names_collapsed <- list(1:(length(new.names)))
   cc <- 1
   for (i in unique(new.names)){
      names_collapsed[cc] <- paste(unlist(new.names[cc]), collapse='')
      cc <- cc + 1
   }
   names(genome)<- names_collapsed
   return(genome)
}


renameGenome(genome, annot)
