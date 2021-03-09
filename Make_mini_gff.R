library(Biostrings)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(parallel)
library(rhdf5)

gff <- 'Schizosaccharomyces_pombe_all_chromosomes.gff3.gz'
annot <- readGFFAsGRanges(gff)
filter_seq <- c("type:CDS", "type:mRNA")

filterGFFByCriteria <- function(gff,criteria)
{
  criteria <- strsplit(unlist(strsplit(criteria,",")),":")
  df_criteria <- as.data.frame(do.call(rbind,criteria))
  criteria_keys <- unique(df_criteria[,1])
  for (i in criteria_keys)
  {
    criteria_conditions <- df_criteria[which(df_criteria[,1] == i),2]
    criteria_conditions_not_na <- criteria_conditions[which(criteria_conditions != "notNA")]
    gff <- gff[mcols(gff)[[i]] %in% criteria_conditions_not_na]
    ## If want to also remove NA values, need to do this using is.na
    if (length(criteria_conditions_not_na) < length(criteria_conditions))
    {
      gff <- gff[!is.na(mcols(gff)[[i]])]
    }
  }
  return(gff)
}

annot <- filterGFFByCriteria(gff = annot,criteria = filter_seq)
unlist_annot <- unlist(GRangesList((annot))) 
annot_df <- as.data.frame(unlist_annot)
small_annot <- annot_df[6:9,]
Fil1_annot <- subset(annot_df, annot_df$Parent %in% c('SPCC1393.08', 'SPCC1393.08.1', 'SPCC1393.08.*'))

result_df <- rbind(small_annot, Fil1_annot)
mini_gff <- GRanges(result_df)

export.gff3(mini_gff, con=file.path('.','mini_Spombe_fil1_and_negative.gff3'))

