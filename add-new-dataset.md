# Setting up a new dataset 

# Table of contents

* [Adding a dataset for a new species](#newspecies)
* [Adding a dataset for an existing species](#existingspecies)


<a name="newspecies"/>


# Adding a dataset for a new species

filler text 


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

