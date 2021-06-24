# Setting up a new dataset 

# Table of contents

* [Setting up a new species](#newspecies)
* [Adding a dataset for an existing species](#existingspecies)


<a name="newspecies"/>


# Setting up a new species

**For genus where annotations are not already available on riboviz**

Some species already having annotations set up in example-datasets, meaning it is quicker and easier to add new datasets for them. The species already present in example-datasets include; *Escherichia*, *Candida*, *Cryptococcus*, *Saccharomyces* and many more. If the species you want to study is not already on example-datasets then it will take more time to set up, but the key steps are outlines below. 

**Contents:**

- [ ] Identify the new species you want to study.
- [ ] Create a genus folder in example-datasets.
- [ ] Download contaminants fasta file.
- [ ] Download or create annotation fasta and gff files.
- [ ] Check annotation files
- [ ] Add your first dataset.

**Identify the new species you want to study**

First step is to pick a species to study. https://www.ncbi.nlm.nih.gov/genome/ contains the genomes of over 60,000 organsims and can provide links to external resources focussing of individual species, for example; Saccharomyces has the https://www.yeastgenome.org/ and schizosaccharomyces pombe has https://www.pombase.org/ listed as externaml resources. 

**Create a genus folder in example-datasets**

Each genus has it's own folder within the relevant Kingdom folder, as outlined in the example-datasets [README file](https://github.com/riboviz/example-datasets#repository-structure-is-loosely-phylogenetic). The genus folder will eventually hold the config.yamls for each dataset. Within the genus folder, an annotation and a contaminants folder need to be created. These will hold the transcriptome gff and fasta, and the contaminant fasta respectively, along with the relevant provinance files.  

**Download contaminants fasta file**

The contaminants file is a fasta file containing the sequences of the rRNA for the organism being studied. rRNA is protected by the ribosome during the digestion step of ribosome profiling, so will be present in the data. These can typically be found and downloaded from a genome database in the .fasta format.

**Download or create annotation fasta and gff files**

Riboviz runs using transcriptome style annotation files rather than genome style annotation files, as this reduces problems later on in the data analysis that would be cause by splicing. If transcriptome fasta and gff files are available to download then great! If not, then these can be created from genome style fasta and gff files using the Rscript 'create_files_for_riboviz.R' available [here](https://github.com/riboviz/riboviz/blob/create_riboviz_style_cds_gff_acope3-278/rscripts/create_files_for_riboviz.R). This code may need to be adapted to fit different research questions and is still in development, but is definitely a good starting point. 

this can be run using the command:

`Rscript --vanilla rscripts/create_files_for_riboviz.R `

and by adding arguments to specify the input fasta and gff, output files, the method of seperating sequences (ie by the gene name or sequence ID) etc.
```
"-i","--input",help="Input DNA sequences. Should contain genome (i.e. sequence of each chromosome) or transcripts. Should be file path",type="character"
"-g","--gff",help="GFF3 file corresponding to input DNA sequences.",type="character"
"--out_dir",help="Output directory for Riboviz-style CDS and GFF files. Will be created if does not exist.",type="character",default="./"
"--out_cds",help="Name of Riboviz-style CDS file",type="character",default="riboviz_cds.fasta"
"--out_gff",help="Name of Riboviz-style GFF3 file",type="character",default="riboviz_gff.gff3"
"-s","--seq_id",help="GFF column to use as the sequence ids",type="character",default="Name"
"-b","--buffer",help="Buffer to use for UTRs",type="integer",default=250
"--h5_file",help="File name for createing H5 file. If not initialized, file will not be created.",type="character",default=NULL
"--codon_data_file",help="File name for codon position .Rdata file. If not initialized, file will not be created.",type="character",default=NULL
"--num_cores",help="Number of cores to use for parallelizable processes.",type="integer",default=1
"--codons_exclude",help="Exclude the first n codons when creating codon_data_file, where n is specified by this argument",default=0
"--remove_trailing",help="Remove trailing info from names to be used for CDS, e.g. remove anything after '_' or '.'",type="character",default="_|\\."
"--filter_seq",help="A comma-separated list of filtering criteria to apply to the GFF3 file, e.g. 'type:CDS,orf_classification:Verified,orf_classification:Uncharacterized'. Use 'notNA' to filter values that are NA, e.g. 'orf_classification:!NA'.",type="character",default="type:CDS"
"--exons_preordered",help="Some GFF3 files have exons pre-ordered such that exon with start codon is listed first. Effects how multi-exon genes will be combined.",action="store_true"

```

Key qualities of transcriptome style fasta and gff files produced by create_ files for riboviz.R:
- All sequences are +ve stranded 
- All coding sequences are flanked by a consistent buffer region
- Each transcript sequence is listed seperately in the fasta file, rather than being the full sequences of chromosomes

**Check annotation files**

**Add your first dataset**


<a name="existingspecies"/>

# Adding a dataset for an existing species 


**For existing species** 

Components: 

- [ ] Identify paper and associated dataset - Link paper
- [ ] Identify the adapter sequence. 
- [ ] Link/describe what the adapter sequence is and where you found it. 
- [ ] Confirm/deny presence of UMIs and barcodes*. 
- [ ] List the UMIs and/or barcodes and link/describe where you found them. Default is FALSE.
- [ ] Identify the ribosome profiling samples from the dataset (some may be RNA-seq)**
- [ ] List the samples and their accession numbers and link to the relevant genome database (e.g. Sequence Read Archive (SRA) or GEO)
- [ ] Using the gathered information, create a config.yaml file***.
- [ ] Prepare to run the data through riboviz****. 

*  *Whether or not UMIs and/or barcodes are present can be hard to pinpoint. If in doubt set the UMIs configuration as FALSE, if this is incorrect this will be reflected in irregular output files*. 

** *The dataset can contain ribosome profiling samples and RNA-seq samples, make sure the ones you select are relevant*. 

*** *It is recommended to use an existing config.yaml file as a guide. Helpful information on parameters and their defaults can be found [here.](https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md)* 

**** *This includes setting up an input directory to hold the sample data in. Details of how to do this can be found [here.](https://github.com/riboviz/riboviz/blob/main/docs/user/run-on-eddie.md#run-a-full-size-example-dataset)*

Testing:

Attempt a run using the riboviz nextflow 
- [ ] Common issues include running the command from the wrong directory, incorrect adapter sequences, UMIs or pathway errors. 
- [ ] Check the logs to see where it failed. 

If the run is successful output files will be generated, check these for:
- [ ] Appropriate lengths – usually 28-32 nt 
- [ ] 3nt periodicity 

If the output files look as expected this can confirm that the run was processed correctly.  

