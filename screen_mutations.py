# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 23:57:19 2026

@author: art_s
"""
# Import libraries
import pandas as pd
from Bio import SeqIO


# Import data
data = {}

for record in SeqIO.parse("mtDNA_seq.fas", "fasta"):
    data[record.id] = str(record.seq)

print(data)

mut_db = pd.read_csv('mtDNA_mutation_database.csv')


# QC filtering
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


# Mutation screening
# Funtional genertation
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


# Filtering only alternate
alt_mut = []
primary = []

for i in results:
    
    if i['Match'] == 'ALT':
        alt_mut.append(i)

        if i.get('Classification') == 'Primary mutations':
            primary.append(i)

# Export analysis results
df_results = pd.DataFrame(results)
df_results.to_csv('output_report.csv', index = False) 

df_alt = pd.DataFrame(alt_mut)
df_alt.to_csv('output_alternate.csv', index = False) 
    










