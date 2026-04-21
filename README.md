# mtDNA-mutation-screening
Python-based mtDNA mutation screening pipeline with QC and clinical annotation, developed as a bioinformatics learning project.
## About the project
Mitochondria are an important cell organelle that functions in generating ATP for cells as a source of chemical energy. This mechanism is involved in several genes in the nucleus and mitochondria. However, the abdominal genes in mitochondria can affect and cause diseases in humans, including Leigh syndrome, Leber hereditary optic neuropathy (LHON), myoclonic epilepsy with ragged red fibers (MERRF), and mitochondrial encephalomyopathy, lactic acidosis, and stroke-like episodes (MELAS). MtDNA diseases result in severe multisystem complications in all ages. In this project, the Python-based pipeline was performed to identify the mutations in mitochondrial DNA (mtDNA) sequences of humans available in the NCBI database. The information on mutations and clinical annotation was from the literature of the Mustafa et al. study.

## Features

## Limitations
* Only supports SNP detection (no insertion/deletion handling)
* Assumes sequences are aligned to rCRS
* Sequences with abnormal length are excluded during QC
## Note
This project was developed as a bioinformatics learning project. The code was written and refined by the author, with assistance from AI tools for code suggestions, debugging, and optimization. Workflow design, data processing, and biological interpretation were performed by the author.
