suppressMessages(library(Rsamtools))
suppressMessages(library(rtracklayer))
suppressMessages(library(rhdf5))
suppressMessages(library(parallel))
suppressMessages(library(optparse))
suppressMessages(library(RcppRoll))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))
suppressMessages(library(here))

## functions to be used in code - replace with source(here::here("rscripts", "read_count_functions.R"))

#Get a data matrix for the gene of interest 
GetGeneDatamatrix <- function(gene, dataset, hd_file){
  rhdf5::h5read(file = hd_file, name = paste0("/", gene, "/", dataset, "/reads/data")) %>%
    return()
}


# read in the gff so positions of start codons and UTRs/buffers can be identified
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

# use the gff to pull out the location of the start codon
GetCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(start) %>%  # pull() pulls out single variable
    min 
}

# use the gff to pull out the location of the stop codon 
GetCDS3end <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(end) %>% 
    max 
}

# Create a datamatrix 
GetGeneDatamatrix5start <- function(gene, dataset, hd_file, 
                                    posn_5start, n_buffer, nnt_gene,posn_3end=Inf) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hd_file)
  
  # if n_buffer bigger than length n_utr5, pad with zeros:
  if (posn_5start > n_buffer) {
    # if posn_5start bigger than n_buffer
    n_left5 <- posn_5start - n_buffer # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # if length n_utr5 less than n_buffer
    n_left5 <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(
      0, 
      nrow = nrow(data_mat_all), 
      ncol = (n_buffer - posn_5start + 1 )
    )
  }
  n_right3 <- posn_5start + nnt_gene - 1 # column to end with (3'end)
  if (n_right3 > posn_3end) {
    
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = n_right3 - posn_3end)
    n_right3 <- posn_3end
  } else{
    zeropad3_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  }
  
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  x<-do.call("cbind",list(zeropad5_mat, data_mat_5start,zeropad3_mat))
  return(x)
}

TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
   # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}
file_url <- "wt.noAT.ribo.4_s.h5"

info <- h5readAttributes(file_url, '/SPCC1393.08.1/D-Sp_2018/reads')

# get the gff as a dataframe - function 



# get the S_pombe gff

gff_df <- readGFFAsDf('S_pombe_full_UTR_or_50nt_buffer.gff3')

# get the cds start 


Fil1_start <- GetCDS5start("SPCC1393.08.1", gff_df)

# get end of sequence 



Fil1_end <- GetCDS3end("SPCC1393.08.1", gff_df)
Fil1_Read_lens <- info$reads_by_len
Fil1_buffer_left <- info$buffer_left
Fil1_buffer_right <- info$buffer_right



Fil1_gene_matrix <- GetGeneDatamatrix5start('SPCC1393.08.1',
                        'D-Sp_2018',
                        file_url,
                        n_buffer = Fil1_buffer_left,
                        nnt_gene = Fil1_end + Fil1_buffer_right,
                        posn_5start = Fil1_start,
                        posn_3end = Fil1_end)

nnt_buffer <- Fil1_buffer_left


Fil1_gene_tidy_start <- TidyDatamatrix(data_mat = info, startpos = -nnt_buffer +1 )
# make tidy data matrix, with all the counts of different read lengths at different positions
# aim: make a tidyDataMatrix with two columns: position and count
# input: Fil1TidyDataMatrix 
# process: for position i, take the sum of the counts at that position across read lengths
#          and return a new row for the new matrix with the position and the counts 
# for each position, take all the counts of reads of different sizes, and combine them into one 
# 'total read' column for each position in the gene
# output: TidyDataMatrix

Fil1_gene_Total_reads_at_position <- tibble(Pos =integer(),
                                            Counts = integer())

for(i in Fil1_gene_tidy_start[1,]$Pos:max(Fil1_gene_tidy_start$Pos)){
  tmp_row <- Fil1_gene_tidy_start %>% filter(Fil1_gene_tidy_start$Pos ==i)
  new_row <- tibble(Pos = i, Counts = sum(tmp_row$Counts))
  Fil1_gene_Total_reads_at_position <- Fil1_gene_Total_reads_at_position %>% bind_rows(new_row)
  }
  

plot(Fil1_gene_Total_reads_at_position$Pos,Fil1_gene_Total_reads_at_position$Counts)

ggplot(Fil1_gene_Total_reads_at_position,aes(x = Pos, y = Counts))+
  geom_col( width = 1, color = 'red')+
  labs(title = 'wt.noAT.ribo.4_s - Fil1', x = 'Position relative to start codon', y = 'Number of reads')
  
 

