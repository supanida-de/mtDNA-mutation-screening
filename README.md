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
```
#### Export QC report
```
df = pd.DataFrame(fail)
df.to_csv('fail_report.csv', index = False)
```
#### Mutation screening
```
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

results = screen(passed, mut_db)
```
#### Filtering only alternate
```
alt_mut = []
patho = []

for i in results:
    
    if i['Match'] == 'ALT':
        alt_mut.append(i)

        if i.get('Classification') == 'Primary mutations':
            patho.append(i)
```
#### Export analysis results
```
df_results = pd.DataFrame(results)
df_results.to_csv('output_report.csv', index = False) 

df_alt = pd.DataFrame(alt_mut)
df_alt.to_csv('output_alternate.csv', index = False) 
```
## Limitations
* Only supports SNP detection (no insertion/deletion handling)
* Assumes sequences are aligned to rCRS
* Sequences with abnormal length are excluded during QC

## Note
This project was developed as a bioinformatics learning project. The code was written and refined by the author, with assistance from AI tools for code suggestions, debugging, and optimization. Workflow design, data processing, and biological interpretation were performed by the author.
