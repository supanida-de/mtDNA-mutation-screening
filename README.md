# mtDNA mutation screening
Python-based mtDNA mutation screening pipeline with QC and clinical annotation, developed as a bioinformatics learning project.

## About the project
Mitochondria are an important cell organelle that functions in generating ATP for cells as a source of chemical energy. This mechanism is mediated by proteins encoded by both nuclear and mitochondrial genes [1-2]. However, the abnormal genes in mitochondria can affect and cause diseases in humans, including Leigh syndrome, Leber hereditary optic neuropathy (LHON) [3], myoclonic epilepsy with ragged red fibers (MERRF) [4], and mitochondrial encephalomyopathy, lactic acidosis, and stroke-like episodes (MELAS) [5]. This project implements a Python-based pipeline to identify mutations in mitochondrial DNA (mtDNA) sequences of humans available in the NCBI database. Mutation information and clinical annotation were curated from the Mustafa et al. study [6].

## Features
* SNP-based mtDNA mutation detection from FASTA sequences  
* Quality control (QC) filtering of input sequences  
* Clinical annotation of mutations (gene, disease, classification)

## Getting Started
### Data
Example input and output datasets are included in this repository.
#### Input
* `mtDNA_seq.fas` — mtDNA sequences in FASTA format
* `mtDNA_mutation_database.csv` — mutation annotation dataset
#### Output
* `output_report.csv` — comprehensive mutation screening results, including both reference and alternate allele matches
* `output_alternate.csv` — filtered results containing only detected variants (alternate allele matches)
* `fail_report.csv` — sequences flagged during quality control (QC)

## How to run
### Dependencies
* Python 3.x
* Pandas
* Biopython
### Installation dependencies
#### Using pip
Install all required packages:
```
# Install Biopython
pip install pandas biopython

# Install pandas
pip install pandas
```
#### Using Anaconda
If you are using Anaconda, `pandas` is typically pre-installed. You only need to install Biopython:
```
pip install biopython
```
### Run the pipeline
#### Import libraries
```
import pandas as pd
from Bio import SeqIO
```
#### Load sequence data
`Bio.SeqIO.parse()` is the core function in the Biopython library used to read sequence data from FASTA files and other file formats. 
```
data = {}

for record in SeqIO.parse("mtDNA_seq.fas", "fasta"):
    data[record.id] = str(record.seq)

print(data)
```
#### Load mutation database
The `pd.read_csv()` function from pandas is used to load mutation data stored in CSV format.
```
mut_db = pd.read_csv('mtDNA_mutation_database.csv')
```
#### QC filtering of sequences
This step performs quality control (QC) by filtering sequences based on their length. Only sequences with the expected mitochondrial genome length (16,569 bp) are retained for downstream analysis, while sequences with abnormal lengths are flagged in the `fail` variable. The fail report is then converted into a DataFrame using `pd.DataFrame()` and exported as a CSV file using `.to_csv()`.
```
passed = {}
fail = []

for name, seq in data.items():
    if len(seq) == 16569:
        passed[name] = str(seq)
    
    elif len(seq) < 16569:
        fail.append({'Sample_ID': name,
                     'Lenght': len(seq),
                     'Reason': 'short_length'})
    else:
        fail.append({'Sample_ID': name,
                     'Lenght': len(seq),
                     'Reason': 'long_length'})

print(fail)

df = pd.DataFrame(fail)
df.to_csv('fail_report.csv', index = False)
```
* `passed`: stores sequences that pass QC and will be used for mutation screening
* `fail`: stores sequences that do not meet the expected length criteria, along with the reason for exclusion
#### Mutation screening
This step performs mutation screening by comparing each nucleotide in the input sequences with a curated mutation database. For each position, the observed base in the sequence is compared against the reference and alternate alleles.
```
# Functional generation
def screen(passed, mut_db):
    results = []
    
    for name, seq in passed.items():
        for index, row in mut_db.iterrows():
            pos = row['Position']
            ref = row['Reference']
            alt = row['Alternate']
            
            seq_index = pos - 1  # convert 1-based → 0-based
            
            if seq_index >= len(seq):
                continue
            
            base = seq[seq_index]
            
            if base == alt:
                results.append({
                    'Sample_ID': name,
                    'Position': pos,
                    'Reference': ref,
                    'Alternate': alt,
                    'Match': "ALT",
                    'Gene': row['Gene'],
                    'Disease': row['Disease'],
                    'Classification': row['Classification']
                })
                
            elif base == ref:
                results.append({
                    'Sample_ID': name,
                    'Position': pos,
                    'Reference': ref,
                    'Alternate': alt,
                    'Match': "REF",
                    'Gene': row['Gene'],
                    'Disease': row['Disease'],
                    'Classification': row['Classification']
                })
    
    return results

# Analysis
results = screen(passed, mut_db)
```
* Iterates through each sequence that passed QC
* Compares each position with the mutation database
* Converts genomic coordinates from 1-based to Python 0-based indexing
* Classifies matches as:
  * ALT: observed base matches the mutation (variant detected)
  * REF: observed base matches the reference (no mutation)
* Stores mutation details along with gene, disease, and classification information
#### Filtering only alternate
This step filters the screening results to retain only detected mutations (alternate allele matches). It also extracts clinically relevant variants classified as primary mutations, which are the most common mutation associated with the disease.
```
alt_mut = []
primary = []

for i in results:
    
    if i['Match'] == 'ALT':
        alt_mut.append(i)

        if i.get('Classification') == 'Primary mutations':
            primary.append(i)
```
* Filters results to include only `ALT` matches (detected variants)
* Stores all detected mutations in `alt_mut`
* Further filters mutations classified as primary into `primary`
* Enables separation of general variants from clinically significant mutations
#### Export analysis results
The results will be generated as CSV files:
* output_report.csv
* output_alternate.csv
```
df_results = pd.DataFrame(results)
df_results.to_csv('output_report.csv', index = False) 

df_alt = pd.DataFrame(alt_mut)
df_alt.to_csv('output_alternate.csv', index = False) 
```

## Limitations
* Only supports SNP detection (no insertion/deletion handling)
* Assumes sequences are aligned to Revised Cambridge Reference Sequence (rCRS) of the human mtDNA
* Sequences with abnormal length are excluded during QC

## Abbreviations
* mtDNA — mitochondrial DNA  
* REF — reference allele  
* ALT — alternate allele (variant)  
* QC — quality control  
### Disease abbreviations
* CPEO — Chronic progressive external ophthalmoplegia
* KSS — Kearns–Sayre syndrome
* LHON — Leber hereditary optic neuropathy
* LS — Leigh syndrome
* MELAS — Mitochondrial myopathy, encephalopathy, lactic acidosis, and stroke-like episodes
* MERRF — Myoclonic epilepsy associated with ragged red fibers
* MIDD — Maternally inherited diabetes and deafness (MIDD)
* NARP — Neurogenic muscle weakness, ataxia, and retinitis pigmentosa

## Author’s note
This project was developed as a bioinformatics learning project. The code was written and refined by the author, with assistance from AI tools for code suggestions, debugging, and optimization. Workflow design, data processing, and biological interpretation were conducted by the author.

## References
[1] Sato, M., and Sato, K. (2013). Maternal inheritance of mitochondrial DNA by diverse mechanisms to eliminate paternal mitochondrial DNA. Biochim Biophys Acta 1833(8), 1979-1984. doi: 10.1016/j.bbamcr.2013.03.010.  
[2] Jonckheere, A.I., Smeitink, J.A., and Rodenburg, R.J. (2012). Mitochondrial ATP synthase: architecture, function and pathology. J Inherit Metab Dis 35(2), 211-225. doi: 10.1007/s10545-011-9382-9.  
[3] Cwerman-Thibault, H., Augustin, S., Ellouze, S., Sahel, J.A., and Corral-Debrinski, M. (2014). Gene therapy for mitochondrial diseases: leber hereditary optic neuropathy as the first candidate for a clinical trial. C R Biol 337(3), 193-206. doi: 10.1016/j.crvi.2013.11.011.  
[4] Mancuso, M., Petrozzi, L., Filosto, M., Nesti, C., Rocchi, A., Choub, A., et al. (2007). MERRF syndrome without ragged-red fibers: the need for molecular diagnosis. Biochem Biophys Res Commun 354(4), 1058-1060. doi: 10.1016/j.bbrc.2007.01.099.  
[5] Wang, Y.X., and Le, W.D. (2015). Progress in diagnosing mitochondrial myopathy, encephalopathy, lactic acidosis, and stroke-like episodes. Chin Med J (Engl) 128(13), 1820-1825. doi: 10.4103/0366-6999.159360.  
[6] Mustafa, M.F., Fakurazi, S., Abdullah, M.A., and Maniam, S. (2020). Pathogenic mitochondria DNA mutations: current detection tools and interventions. Genes (Basel) 11(2). doi: 10.3390/genes11020192.
