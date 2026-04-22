# mtDNA mutation screening
Python-based mtDNA mutation screening pipeline with QC and clinical annotation, developed as a bioinformatics learning project.

## About the project
Mitochondria are an important cell organelle that functions in generating ATP for cells as a source of chemical energy. This mechanism is involved in several genes in the nucleus and mitochondria. However, the abnormal genes in mitochondria can affect and cause diseases in humans, including Leigh syndrome, Leber hereditary optic neuropathy (LHON), myoclonic epilepsy with ragged red fibers (MERRF), and mitochondrial encephalomyopathy, lactic acidosis, and stroke-like episodes (MELAS). This project implements a Python-based pipeline to identify mutations in mitochondrial DNA (mtDNA) sequences of humans available in the NCBI database. Mutation information and clinical annotation were curated from the Mustafa et al. study.

## Features
* SNP-based mtDNA mutation detection from FASTA sequences  
* Quality control (QC) filtering of input sequences  
* Clinical annotation of mutations (gene, disease, classification)

## Getting Started
### Data
Example input and output datasets are included in this repository.
#### Input
* mtDNA_seq.fas
* mtDNA_mutation_database.csv
#### Output
* fail_report.csv
* output_alternate.csv
* output_report.csv 

## Limitations
* Only supports SNP detection (no insertion/deletion handling)
* Assumes sequences are aligned to rCRS
* Sequences with abnormal length are excluded during QC

## Note
This project was developed as a bioinformatics learning project. The code was written and refined by the author, with assistance from AI tools for code suggestions, debugging, and optimization. Workflow design, data processing, and biological interpretation were performed by the author.
