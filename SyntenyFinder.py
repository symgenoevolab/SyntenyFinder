#!/usr/bin/env python
# coding: utf-8

# In[5]:


print('Importing dependencies')
import subprocess 
import os
import re
import zipfile
import pandas as pd
import argparse
import concurrent.futures
from Bio import SeqIO 
from io import StringIO
from dependencies.Synteny_functions import Synteny, NCBI_Synteny, run_NCBI_Synteny


# In[2]:


def SyntenyFinder(): 
 

    parser = argparse.ArgumentParser(description = 'Run synteny analysis')
    parser.add_argument('--accessions', 
                        type = str, 
                        help = 'Comma separated list of accessions for NCBI annotated genome assemblies') 
    parser.add_argument('--orthofinder', 
                        type = str, 
                        help = 'Path to executable orthofinder if not installed') 
    parser.add_argument('--threads', '-t', 
                        type = int, 
                        help = 'Number of threads to use when running orthofinder') 
    parser.add_argument('--run_name', '-n', 
                        type = str, 
                        help = 'Name of this run iteration (will become directory name)') 
    parser.add_argument('--directory', '-d', 
                        type = str, 
                        help = 'Root directory for analysis folder to be created') 
    parser.add_argument('--algs', 
                        type = str, 
                        help = 'Species code or accession for dataset used to trace ancestral lineage groups') 

    args = parser.parse_args() 
    
    if args.run_name is not None: 
        print(f'Run name: {args.run_name}') 
    else: 
        raise Exception('A name must be provided for this run iteration')

    if args.accessions is not None: 
        accessions = args.accessions.split(',')
        print(f'Accessions: {accessions}') 
    else: 
        raise Exception('List of accesssions for NCBI annotated genome assemblies must be provided')

    if args.algs is not None: 
        print(f'Tracing lineage groups from chromosomes of species code {args.algs}')
        algs = args.algs
    else: 
        print(f'No alg species provided. Tracing groups from chromosomes of first species listed: {accessions[0]}') 
        algs = accessions[0]
        
    if args.directory is not None: 
        print(f'Root directory: {args.directory}')
        directory = args.directory
    else: 
        print(f'No root directory passed. Will run rooted in current location.') 
        directory = './'
        
    if args.orthofinder is not None: 
        print(f'Orthofinder location: {args.orthofinder}')
        orthofinder = args.orthofinder
    else: 
        print(f'No file location specified for executable orthofinder. Will assume installed.')
        orthofinder = 'orthofinder'

    if args.threads is not None: 
        print(f'Threads: {args.threads}')
    else: 
        print(f'Number of threads unspecified. Will run with OrthoFinder default.') 

    
    synteny = run_NCBI_Synteny(accessions = accessions, 
                               run_name = args.run_name, 
                               algs = algs, 
                               root_directory = directory, 
                               orthofinder_path = orthofinder,
                               threads = args.threads)
                               
if __name__ == "__main__": 
    SyntenyFinder() 


# # In[6]:


# get_ipython().run_line_magic('run', 'dependencies/Synteny_functions.ipynb')
# synteny = run_NCBI_Synteny(accessions = ['GCF_000003815.2','GCF_902652985.1'] , 
#                            run_name = 'Bfl_Pma_sample', 
#                            root_directory = './', algs = 'first', 
#                            threads = 100)






