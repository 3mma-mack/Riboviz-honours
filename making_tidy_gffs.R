
library(tidyverse)
library(rtracklayer)

genome <- readDNAStringSet('S_pombe_full_UTR.fasta',format = "fasta")
annot <- readGFF('Alex_files/Schizosaccharomyces_pombe_CDS_w_250utrs.gff3')

test <- tibble(annot)

three_UTR <- readDNAStringSet('3UTR.fa.gz')
CDS <- readDNAStringSet('cds.fa.gz')
five_UTR <- readDNAStringSet('5UTR.fa.gz')

names(five_UTR) <- paste(names(five_UTR), '.1', sep = '')
names(CDS) <- paste(names(CDS), '.1', sep = '')
names(three_UTR) <- paste(names(three_UTR), '.1', sep = '')

# seperate out the UTRs and CDS to allow processing 

gff_fiveutr <- filter(test, test$type == 'UTR5')
gff_cds <- filter(test, test$type == 'CDS')
gff_threeutr <- filter(test, test$type == 'UTR3')

i <- 'SPAC1002.02.1' 

# Edit the ranges of the 5UTR gff to be 1-width of listed UTR for each gene, if present

gff_fiveutr_edit <- tibble(seqid = factor(),
                           source = factor(),
                           type = factor(),
                           start = integer(),
                           end = integer(),
                           score = double(),
                           strand = character(),
                           phase = integer(),
                           Name = character())

for(i in gff_fiveutr$Name){
  if(i %in% names(CDS)){
    if(i %in% names(five_UTR)){
      new_row <- filter(gff_fiveutr, gff_fiveutr$Name == i) %>%
        mutate(start = 1, end = width(five_UTR[names(five_UTR)==i])) 
      gff_fiveutr_edit <- gff_fiveutr_edit %>% bind_rows(new_row)
      }else{
         new_row <- filter(gff_fiveutr, gff_fiveutr$Name == i) %>%
           mutate(start = 0,end = 0) 
          gff_fiveutr_edit <- gff_fiveutr_edit %>% bind_rows(new_row)       
    }
  }
}


# repeat for CDS, checking if a 5' UTR is present 
gff_cds_edit <- tibble(seqid = factor(),
                       source = factor(),
                       type = factor(),
                       start = integer(),
                       end = integer(),
                       score = double(),
                       strand = character(),
                       phase = integer(),
                       Name = character())

for(i in gff_cds$Name){
  if(i %in% names(CDS)){
    if(i %in% names(five_UTR)){
      new_row <- filter(gff_cds, gff_cds$Name == i) %>%
        mutate(start = width(five_UTR[names(five_UTR)==i]) +1,
               end = width(CDS[names(CDS)==i]) + width(five_UTR[names(five_UTR)==i]) +1)
      gff_cds_edit <- gff_cds_edit %>% bind_rows(new_row)
      
    }else{
      new_row <- filter(gff_cds, gff_cds$Name == i) %>%
        mutate(start = 1,
               end = width(CDS[names(CDS)==i]))
      gff_cds_edit <- gff_cds_edit %>% bind_rows(new_row)
    }
  }
}

gff_threeutr_edit <- tibble(seqid = factor(),
                       source = factor(),
                       type = factor(),
                       start = integer(),
                       end = integer(),
                       score = double(),
                       strand = character(),
                       phase = integer(),
                       Name = character())

for(i in gff_threeutr$Name){
  if(i %in% names(CDS)){
    if(i %in% names(three_UTR)){
      new_row <- filter(gff_threeutr, gff_threeutr$Name == i) %>%
        mutate(start = (filter(gff_cds_edit, gff_cds_edit$Name == i)$end +1),
             end = (filter(gff_cds_edit, gff_cds_edit$Name == i)$end + 
               width(three_UTR[names(three_UTR)==i]) +1 ))
      gff_threeutr_edit <- gff_threeutr_edit %>% bind_rows(new_row)
  
      }else{
        new_row <- filter(gff_threeutr, gff_threeutr$Name == i) %>%
          mutate(start = (filter(gff_cds_edit, gff_cds_edit$Name == i)$end),
             end = (filter(gff_cds_edit, gff_cds_edit$Name == i)$end))
        gff_threeutr_edit <- gff_threeutr_edit %>% bind_rows(new_row)
      }
   }
}


total_gff <- tibble(seqid = factor(),
                    source = factor(),
                    type = factor(),
                    start = integer(),
                    end = integer(),
                    score = double(),
                    strand = character(),
                    phase = integer(),
                    Name = character())

for(i in 1:nrow(gff_cds_edit)){
  total_gff<- total_gff %>% rbind(gff_fiveutr_edit[i,])
  total_gff<- total_gff %>% rbind(gff_cds_edit[i,])
  total_gff<- total_gff %>% rbind(gff_threeutr_edit[i,])
}

export.gff3(total_gff, con=file.path('.','S_pombe_tiddle_full_UTR.gff3'))
