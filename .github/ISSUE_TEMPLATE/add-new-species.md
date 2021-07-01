---
name: Add new species
about: 'guide for adding a new species to example-datasets. '
title: ''
labels: ''
assignees: ''

---

- [ ] Identify the new species you want to study.
- [ ] Create a new folder in example-datasets. Outlines on folder structure can be found [here](https://github.com/riboviz/example-datasets#repository-structure-is-loosely-phylogenetic) 
- [ ] Create Annotation and Contaminants folders
- [ ] Download contaminants fasta file. This should contain the rRNA sequences for your species, and (optionally) the tRNA sequences.
- [ ] Add contaminants fasta file to Contaminants folder
- [ ] If available, download transcriptome style annotation fasta and gff files. 
- [ ] If only genome annotation files are available, transcriptome annotation files can made using [create_files_for_riboviz.R](https://github.com/riboviz/riboviz/blob/create_riboviz_style_cds_gff_acope3-278/rscripts/create_files_for_riboviz.R) 
- [ ] Check annotation files using check_fasta_gff. Documentation on use and outputs can be found [here](https://github.com/riboviz/riboviz/blob/main/docs/user/check-fasta-gff.md)
- [ ] Add annotation files to Annotation folder
- [ ] Create a provenance.txt file.
- [ ] Add your first dataset.
