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

- Create a new branch in example datasets to work in.
- Within the branch, create a folder for your genus.
- Inside the genus folder, create an annotation and a contaminants folder. 

**Download contaminants fasta file**

The contaminants file is a fasta file containing the sequences of the rRNA for the organism being studied. rRNA is protected by the ribosome during the digestion step of ribosome profiling, so will be present in the data. These can typically be found and downloaded from a genome database in the .fasta format.

This can then be added to the contaminants folder, along with a provenance.txt file describing the origins and how it was created. 

- Find and download the rRNA sequences for your organisim in the .fasta format
- Add the rRNA.fasta file to the contaminants file 
- Add a provenance.txt file

**Download or create annotation fasta and gff files**

Riboviz runs using transcriptome style annotation files rather than genome style annotation files, as this reduces problems later on in the data analysis that would be cause by splicing. If transcriptome fasta and gff files are available to download then great! If not, then these can be created from genome style fasta and gff files using the Rscript 'create_files_for_riboviz.R' available [here](https://github.com/riboviz/riboviz/blob/create_riboviz_style_cds_gff_acope3-278/rscripts/create_files_for_riboviz.R). This code may need to be adapted to fit different research questions and is still in development, but is definitely a good starting point. 

This can be run using the command:

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
These parameters allow the code to be adapted depending on the format of the GFF and fasta files, which may vary based on how they were created.

Key qualities of transcriptome style fasta and gff files produced by create_ files for riboviz.R:
- All sequences are +ve stranded 
- All coding sequences are flanked by a consistent buffer region
- Each transcript sequence is listed seperately in the fasta file, rather than being the full sequences of chromosomes
- For each transcript, the GFF will contain 3 features; the upstream buffer, the CDS and the downstream buffer.

Once created, these annotation files can be uploaded to the annotation folder. 

- If available, download transcriptome gff and fasta files for the species of interest and add to annnotation folder.
- If unavailable, download genome gff and fasta files and run through create_files_for_riboviz.R
- Add new files to annotation folder

**Check annotation files**

The GFF file acts as a guide, and tells the pipeline where to find the start and stop codons for each transcript. However it is important that these locations do actually match the location of the start and stop codon in the provided fasta file. Files can be checked using the riboviz tool check_fasta_gff, which looks for start and stop codons in the fasta file using the positions provided by the GFF and provides details on the number of missing and unexpected features. Details about how to use check_fasta_gff can be found [here](https://github.com/riboviz/riboviz/blob/main/docs/user/check-fasta-gff.md). It is likely that there will be a few issues present due to alternative start and stop codons, and the presence of psuedogenes if not removed.

If there are issues present, then using a genome viewer such as SnapGene can help highlight the cause.

If the gff and fasta files show few issues, the files are ready to test with a dataset, which you can add following the instructions below. Adding a provenance.txt file here describing how the files were made and the sources of the data is advisable. 

- Run GFF and fasta files through check_fasta_gff.
- Check outputs for unexpected issues.
- Investigate issues if present. 
- Add provenance.txt file.

<a name="existingspecies"/>

# Adding a dataset for an existing species 

If the species you would like to study is already available in example-datasets, the necessary annotation and contamination files are already available. You will need to ensure that [riboviz is set up](https://github.com/riboviz/riboviz/blob/main/docs/user/install.md) and that the [vignette runs](https://github.com/riboviz/riboviz/blob/main/docs/user/run-on-eddie.md#run-a-vignette-of-the-riboviz-workflow) before adding another dataset.

This documentation provides the user with some useful documentation which supplements the [documentation](https://github.com/riboviz/riboviz/blob/main/docs/user/run-on-eddie.md#run-a-full-size-example-dataset) available for  new datasets for species which have already been set up.

- [ ] Identify scientific paper and associated dataset - link the paper
- [ ] Identify the ribosome profiling samples from the dataset (some may be RNA-seq) - link database
- [ ] Identify the adapter sequence - link/describe  
- [ ] Confirm/deny presence of UMIs and barcodes  
- [ ] Using the gathered information, create a config.yaml file.
- [ ] Prepare to run the data through riboviz. 

**Components:** 

The scientific paper 
- Identify the scientific paper which is associated with the dataset of interest. The paper will contain a reference to the associated sequencing data, usually in the form of which database the sequencing data has been deposited in and the series **accession number**. 
- The methods section of the paper will either refer to a protocol or directly list the **adapter sequence** used as well as **UMIs** and/or **barcodes** which may have been used. This information is needed for the config file. 

The ribosome profiling samples 
- To find the fastq files for the samples associated with the dataset the series accession number can be searched for in the relevant database.
- Identify the ribosome profiling samples as some of them may be RNA-seq and should not be included. ([Instructions for downloading fastq files](https://github.com/riboviz/riboviz/blob/main/docs/user/run-on-eddie.md#download-fastq-data-files-from-the-short-read-archive-sra-initial-setup))

The config.yaml file 
- Documentation for [configuring the config.yaml](https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md) file is available in addition to a wide variety of examples in [example-datasets](https://github.com/riboviz/example-datasets). It is recommended to use an existing config.yaml file as a guide.
- The majority of the parameters have default settings that are species-specific and can be modelled off of previous config files available in example-datasets. 
- Dataset-specific parameters include the adapter sequence, accession codes for the samples and the potential presence of UMIs and/or barcodes. Whether or not UMIs and/or barcodes are present can be hard to pinpoint. If in doubt set the UMIs configuration as FALSE, if this is incorrect this will be reflected in irregular output files. 
- Once your config.file has been written it is useful for fellow users to be able to access your work. To do this create a branch in example-datasets with a relevant name (e.g. W-Sc_2016_63 for Weinberg et al 2016 Saccharomyces cerevisiae dataset from issue ticket 63) 
- From this branch add the config file and commit the config file with a useful commit message so that other users can easily access your work via this branch in example datasets. It is useful to link/name this branch in your issue ticket so that it is easy to find. 

**Testing:**

Attempt a run using the riboviz nextflow 
- [ ] Common issues include running the command from the wrong directory, incorrect adapter sequences, UMIs or pathway errors. 
- [ ] Check the logs to see where it failed. 

If the run is successful output files will be generated, check these for:
- [ ] Appropriate lengths – usually 28-32 nt 
- [ ] 3nt periodicity 

If the output files look as expected this can confirm that the run was processed correctly.  

