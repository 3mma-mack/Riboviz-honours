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