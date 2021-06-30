---
name: Add new example dataset
about: 'Guide for adding a new dataset to example-datasets. '
title: ''
labels: ''
assignees: ''

---

When adding a new example dataset providing links and listing adapters will make the process easier for others to follow. Key steps are listed below.

- [ ] Identify paper or data source - list and link
- [ ] Identify the ribosome profiling samples from the dataset (some may be RNA-seq) - link database
- [ ] Identify adapter sequence - provide sequence
- [ ] Confirm or deny presence of UMIs and barcodes if used - describe if present
- [ ] If UMIs are present, create UMI regular expression. Examples are described [here](https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md#examples)
- [ ] Using information gathered, create config file. Information on parameters can be found [here](https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-config.md#configuration-parameters) 
- [ ] Download sample data 
- [ ] Set up riboviz, following [documentation](https://github.com/riboviz/riboviz/blob/main/docs/user/prep-riboviz-run-nextflow.md). If using Eddie, documentation is also available [here.](https://github.com/riboviz/riboviz/blob/main/docs/user/run-on-eddie.md)
- [ ] (optional) Create downsampled data and fast test run on that
- [ ] Test run of full sized dataset
- [ ] Look at results - check for 3nt periodicity in coding regions, most common read lengths being 28-32 nt, and clear start and stop profiles
- [ ] Troubleshoot if necessary 
- [ ] Update genus-level README.md and provenance section of config file
