#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Isabel Jiah-Yih Liao
# January 2024 
# Synteny helper functions

# This file contains the helper functions needed to run Synteny. 


# In[ ]:


import subprocess 
import os
import re
import zipfile
import pandas as pd
import argparse
import concurrent.futures
from Bio import SeqIO 
from io import StringIO


# In[7]:


class Synteny: 
    def __init__(self, root_directory, run_name, species_codes, 
                 proteome_ext = '.faa', features_ext = '.gff', genome_ext = '.fna'): 
        self.run_name = run_name 
        self.dirs = {'root' : root_directory} 
        self.metadata = pd.DataFrame(index = species_codes)
        self.species_data = {species:dict() for species in species_codes}
        self.check(proteome_ext, annotation_ext, genome_ext)
        self.build()
        try: 
            self.incorporate_orthofinder()
            print('OrthoFinder results found at ' + self.dirs['orthofinder_results'])
        except: 
            print('No orthofinder results found.')

    ###########
    ### Initialising functions 
    ###########
    # Check to see if all the necessary files are correctly named and can be found. 
    # Returns a dictionary dirs containing the appropriate file paths containing the input data
    def check(self, proteome_ext): 
        self.add_dir('root', 'input_data') 
        for type, extension in [('proteomes', proteome_ext), 
                                ('genomes', genome_ext), 
                                ('features', features_ext)]: 
            self.add_dir('input_data', type)

            # print(f"Checking files in {self.dirs[type]}...")
            contents = os.listdir(self.dirs[type])
            contents = [filename for filename in contents if filename.endswith(extension)]
            for species in self.metadata.index: 
                files = [filename for filename in contents if filename.startswith(species)]
                if len(files) != 1: 
                    print(f"Incorrect number of {type} files found for {species} \n{files}")
                else: 
                    identifier = f"{species}_{type}"
                    self.dirs[identifier] = os.path.join(self.dirs[type], files[0])
        print('File check complete') 
        return

    # Builds neccesary directories to run
    def build(self): 
        
        self.make_dir('root', self.run_name)
        self.make_dir(self.run_name, 'run_files')
        self.make_dir(self.run_name, 'output')
        self.make_dir('run_files', 'run_proteomes') 
        print('Directories built for run...') 
       
        return 
        
    ###########
    ### Preparing proteomes
    ###########    
    # Read proteome files and save the SeqRecords into a list. 
    # Return a dictionary sending each species code to its respective SeqRecords list. 
    def read_proteomes(self): 
        print('Reading proteomes...') 
        for species in self.metadata.index: 
            filepath = self.dirs[f'{species}_proteomes']
            records = []
            with open(filepath, 'r') as handle: 
                for record in SeqIO.parse(handle, 'fasta'): 
                    records.append(record)
            self.species_data[species]['proteomes'] = records
        return self.species_data[species]['proteomes'][:3]

    # For algs provided through adapted gene names in the proteome, parse the 
    # algs and simplify the gene names prior to running orthofinder. 
    def proteome_id_trim(self, species, position, delimiter = '_'): 
        index = position - 1 
        
        species_proteome = self.species_data[species]['proteomes']
        gene_ids = {record.id : record.id.split(delimiter)[index] for record in species_proteome}
       
        for gene in self.species_data[species]['proteomes']: 
            gene.id = gene_ids[gene.id]
            
        return self.species_data[species]['proteomes'][:3]

        
    # Adds the three letter species code before each gene to make orthofinder output
    # easer to understand. 
    # If no list of species codes is given, all proteomes are modified. 
    def proteome_add_species(self, species_list = None): 
        if species_list is None: 
            species_list = self.metadata.index
        try: 
            for species in species_list: 
                for gene in self.species_data[species]['proteomes']: 
                    gene.id = species + '|' + gene.id
                    gene.description = ''
                print('Proteome gene names modified for ' + species)
            return self.species_data[species]['proteomes'][:3]
        except: 
            print('Error occured in adding species codes to proteome.')
            print('Check that read_proteomes() has been run on all species') 
            return 

    # Writes proteomes selected for use to run_proteomes folder in run_files
    def write_proteomes(self): 
        try: 
            for species in self.metadata.index: 
                new_proteome_file = os.path.join(self.dirs['run_proteomes'], f"{species}.faa")
                self.dirs[f"{species}_run_proteomes"] = new_proteome_file 
                with open(new_proteome_file, 'w') as output_handle: 
                    SeqIO.write(self.species_data[species]['proteomes'], output_handle, 'fasta')
                print(f'{species} proteomes successfully written to {new_proteome_file}.')
            return self.species_data[species]['proteomes'][:3]
        except: 
            print('Error occured in writing proteomes. Check that read_proteomes() has been run.') 
            return
            
    # Function to primary_transcript.py as specified in the OrthoFinder tutorial?

    ###########
    ### Orthofinder
    ###########
    # Running orthofinder using nohup 
    def run_orthofinder(self, orthofinder_path = 'orthofinder', threads = None):
        self.dirs['orthofinder_executable'] = orthofinder_path 
        orthofinder_output = os.path.join(self.dirs['run_files'], 'orthofinder_output') 
        if os.path.exists(orthofinder_output): 
            print(f"Previous OrthoFinder output directory detected at {orthofinder_output}." \
                  f"Continuing run on existing data."\
                  f"If you would like to rerun OrthoFinder, please remove the specified directory.") 
            return
        else: 
            command = ['nohup', 
                       orthofinder_path, 
                       '-f', self.dirs['run_proteomes'], 
                       '-o', orthofinder_output]
            if threads is not None: 
                command = command + ['-t', str(threads)]
            print(command)
            
            result = subprocess.Popen(command, stdout=subprocess.PIPE)
            for line in result.stdout: 
                print(line.decode(), end = '')
            result.wait()
            return 

    # Incorporating Orthofinder results
    def incorporate_orthofinder(self): 
        self.add_dir('run_files', 'orthofinder_output') 
        results_dir = os.listdir(self.dirs['orthofinder_output'])[-1]
        self.dirs['orthofinder_results'] = os.path.join(self.dirs['orthofinder_output'], results_dir)
        self.add_dir('orthofinder_results', 'Orthogroups')
        self.dirs['orthogroups'] = os.path.join(self.dirs['Orthogroups'], 'Orthogroups.tsv')
        self.orthogroups = pd.read_csv(self.dirs['orthogroups'], sep = '\t') 
        self.dirs['sco'] = os.path.join(self.dirs['Orthogroups'], 'Orthogroups_SingleCopyOrthologues.txt')
        return

    # Return a table of single copy orthologues for a subset of species
    # While it is possible to use the list of single copy orthologues generated by orthofinder, 
    # this method allows for more flexibility and will result in a larger number of shared genes 
    # if run with only a subset of the species used for orthofinder. 
    def single_copy_orthologues(self, species_list = None): 
        try: 
            if species_list is None: 
                species_list = self.metadata.index.tolist()
            orthogroup_subset = (self.orthogroups.copy()
                                       [['Orthogroup', ] + species_list]).dropna()
            multiple_copies = orthogroup_subset.map(lambda x: ',' in x) 
            self.single_copy_orthogroups  = orthogroup_subset[~multiple_copies.any(axis = 1)] 

        except Exception as e: 
            raise Exception("Species names not found in OrthoFinder output. " \
                            "Please ensure a unique run name is specified for each group of species used "\
                            "such that OrthoFinder can be rerun on all included species.")
            return e
        
        return self.single_copy_orthogroups

    ###########
    ### Building the karyotype file from genome
    ###########   
    # Add karyotype information to metadata, then call build_karyotype. 
    def add_karyotype(self, chromosomes_per_species): 
         self.metadata['Chromosomes'] = self.metadata.index.map(chromosomes_per_species) 
         self.build_karyotype()
    
    def build_karyotype(self): 
        print('Reading karyotype files: this may take a moment...') 
        for species in self.metadata.index: 
            species, karyotype_df = self.species_karyotype(species)
            self.species_data[species]['karyotype'] = karyotype_df
        return

    def species_karyotype(self, species): 
        species_karyotype = self.metadata.loc[species]['Chromosomes'] 
        print(f"Reading {self.metadata.loc[species]['Organism_Name']} genome file...")
        # Using SeqIO to read the genome file, get the scaffold id and length for each entry. 
        genome_path = self.dirs[f'{species}_genomes']
        genome_records = list(SeqIO.parse(genome_path, 'fasta'))
        scaffolds = [{'Scaffold': record.id, 'Length': len(record.seq)} 
                     for record in genome_records] 

        # Sort the scaffolds from longest to shortest and get the top scaffolds
        scaffolds = (pd.DataFrame(scaffolds).sort_values('Length', ascending = False)
                                            .reset_index(drop = True))
        chromosomes = scaffolds.head(species_karyotype) 
        chromosomes.columns = ['Chromosome', 'Length'] 
        return species, chromosomes
        
    def clean_karyotype(self, columns, labels = None): 
        for species in self.metadata.index: 
            karyotype = self.species_data[species]['karyotype']
            cleaned_karyotype = pd.DataFrame(index = range(karyotype.shape[0]))
            for header in columns: 
                if header in karyotype.columns: 
                    cleaned_karyotype[header] = karyotype[header]
                elif header == 'SPECIES': 
                    cleaned_karyotype[header] = species
                else: 
                    cleaned_karyotype[header] = header
            if labels is not None: 
                cleaned_karyotype.columns = labels 
            self.species_data[species]['cleaned_karyotype'] = cleaned_karyotype 
        return cleaned_karyotype

    def write_karyotype(self, index = False): 
        try: 
            for species in self.metadata.index: 
                file_path = os.path.join(self.dirs['output'], f"{species}_karyotype.txt")
                self.species_data[species]['cleaned_karyotype'].to_csv(file_path, 
                                                                       sep = '\t', 
                                                                       index = index) 
            print(f"Karyotype files successfully written to {self.dirs['output']}")
        except: 
            raise Exception('Error writing karyotype dataframes to file') 
                
    ###########
    ### Incorporating gene positions from GTF
    ###########
    # Query to keep only rows with certain key words 
    def gtf_to_dataframe(self, species, annotation_type = None, feature = None): 
        gtf_filepath = self.dirs[f'{species}_features']
        columns = ["sequence", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        gtf_dataframe = pd.read_csv(gtf_filepath, sep = '\t', 
                                    comment = '#', 
                                    header = None, 
                                    names = columns, 
                                    dtype={'start': int, 'end': int})
    
        # Keep only rows of a given feature 
        if feature is not None: 
            gtf_dataframe = gtf_dataframe[gtf_dataframe.feature == feature]
        
        # Extract the specified annotation type from the 'attribute' column
        if annotation_type is None: 
            gtf_dataframe['annotation'] = gtf_dataframe['attribute'] 
        else: 
            gtf_dataframe['annotation'] = (gtf_dataframe['attribute'].str
                                          .extract(f'{annotation_type}[ =]"?([^";]+)"?', expand=False))
            gtf_dataframe = gtf_dataframe.dropna() 
        gtf_dataframe = gtf_dataframe.drop_duplicates(subset = 'annotation', keep = 'first') 
        self.species_data[species]['features'] = gtf_dataframe
        return gtf_dataframe

    # Add a species code with a pipe to the beginning of the chosen annotation names to match proteome gene name format 
    def gtf_add_species(self, species_list = None): 
        if species_list is None: 
            species_list = self.metadata.index
        for species in species_list: 
            all_files_read = True
            if 'features' not in self.species_data[species]: 
                print(f'gtf file for {species} not yet read.') 
                all_files_read = False 
        if all_files_read == True: 
            for species in self.metadata.index: 
                new_annotations = species + '|' + self.species_data[species]['features']['annotation']
                self.species_data[species]['features']['annotation'] = new_annotations
        return new_annotations

    # Merges sco genes with the associated GTF information 
    def merge_gtf(self): 
        for species in self.metadata.index: 
            features = self.species_data[species]['features'].copy()
            orthogenes = self.single_copy_orthogroups[[species]].copy()
            orthogene_coords = pd.merge(orthogenes, features, 
                               left_on = species, right_on = 'annotation', 
                               how = 'left')
            self.species_data[species]['orthogene_coords'] = orthogene_coords
        return orthogene_coords

            
    ###########
    ### Incorporating lineage groups 
    ###########
    # After merging the dataframes and creating the chromosome files, trace the genes based on 
    # their chromosome in a given species
    def trace_chromosomes(self, species): 
        groups = self.species_data[species]['cleaned_karyotype']['Chr'].tolist()
        # print(groups)
        group_ids = [chr(ord('A') + i) for i in range(len(groups))]
        genes_to_groups = dict(zip(groups, group_ids))
        genes = self.species_data[species]['orthogene_coords']['sequence']
        gene_groups = [genes_to_groups[gene] if gene in genes_to_groups.keys() else gene for gene in genes ]
        self.algs = gene_groups 
        return 
        
    # For algs provided through adapted gene names in the proteome, parse the 
    # algs and simplify the gene names prior to running orthofinder. 
    def parse_algs(self, species, delimiter = '_', columns = ['index', 'alg', 'gene_id'], 
                                                               mapping_add_species = True): 
        try: 
            alg_index = columns.index('alg') 
        except: 
            raise Exception('alg column not specified') 
        try: 
            gene_id_index = columns.index('gene_id')
        except: 
            raise Exception('gene_id column not specified') 
        
        filepath = self.dirs[f'{species}_proteomes']
        proteome_headers = []
        with open(filepath, 'r') as handle: 
            for record in SeqIO.parse(handle, 'fasta'): 
                proteome_headers.append(record.id)
        gene_headers = pd.DataFrame([record.split(delimiter) for record in proteome_headers]) 
        gene_headers.columns = columns
        if mapping_add_species: 
            gene_headers['gene_id'] = species + '|' + gene_headers['gene_id']
        self.alg_mapping = dict(zip(gene_headers['gene_id'], gene_headers['alg']))
        return gene_headers
        
    def trace_algs(self, species): 
        if not hasattr(self, 'alg_mapping'): 
            raise Exception('No alg_mapping found. Please ensure parse_alg() has been run.')
        genes = self.species_data[species]['orthogene_coords']['annotation'] 
        algs = [self.alg_mapping[gene] if gene in self.alg_mapping.keys() else gene for gene in genes ]
        self.algs = algs 
        return

    ###########
    ### Format and write coordinate files
    ###########  
    # Reformat the merged coordinates dataframe. 
    def clean_coords(self, columns, labels = None): 
        for species in self.metadata.index: 
            coords = self.species_data[species]['orthogene_coords'] 
            cleaned_coords = pd.DataFrame(index = range(coords.shape[0]))
            for header in columns: 
                if header == "ALG": 
                    cleaned_coords[header] = self.algs
                elif header in coords.columns: 
                    cleaned_coords[header] = coords[header]
                else: 
                    cleaned_coords[header] = header
            self.species_data[species]['cleaned_coords'] = cleaned_coords
            if labels is not None: 
                cleaned_coords.columns = labels
        return cleaned_coords 

    # Write the cleaned coordinates to a csv file 
    def write_coords(self): 
        for species in self.metadata.index: 
            file_path = os.path.join(self.dirs['output'], f"{species}_coordinates.tsv") 
            self.species_data[species]['cleaned_coords'].to_csv(file_path, 
                                                                sep = '\t', 
                                                                header = False)
        print(f"Coordinate files successfully written to {self.dirs['output']}")
   
    ###########
    ### Compatibility helpers
    ###########
    def dotplot(self, x_species, y_species): 
        x_karyotype = self.species_data[x_species]['karyotype']
        y_karyotype = self.species_data[y_species]['karyotype']
        x_coords = self.species_data[x_species]['orthogene_coords']
        y_coords = self.species_data[y_species]['orthogene_coords']
        test = pd.merge(x_coords, y_coords, on = 'annotation')
        return test
    
    ###########
    ### Compatibility helpers
    ###########
    # Given a species code and a substring which indicates a suffix, truncate the gene names in 
    # the single_copy_orthogroups dataframe to omit the suffix. 
    def truncate_sco(self, species, suffix): 
        self.single_copy_orthogroups = (self.truncate(
                                                 self.single_copy_orthogroups, 
                                                      species, suffix))
        return self.single_copy_orthogroups
        
    # Given a species code and a substring which indicates a suffix, truncate the annotations in 
    # the features dataframe to omit the suffix. 
    def truncate_features(self, species, suffix):      
        self.species_data[species]['features'] = (self.truncate(
                                                                 self.species_data[species]['features'], 
                                                                     'annotation', suffix)) 
        return self.species_data[species]['features']

        
    # Save the object to a file so that the analysis can be continuted from a later point
    def save_synteny(self): 
        with open(self.dirs['save_synteny'], 'wb') as file: 
            pickle.dump(self, file)

    ###########
    ### Helper functions 
    ###########
    # Add a directory to dirs, assuming that the parent is already in dirs
    def add_dir(self, parent, add): 
        new_path = os.path.join(self.dirs[parent], add)
        if not os.path.exists(new_path): 
            raise Exception(f"Expected directory at {new_path} was not found.") 
            
        else: 
            self.dirs[add] = new_path
        return new_path

    # Create a new directory and add it to dirs
    def make_dir(self, parent, make): 
        new_dir = os.path.join(self.dirs[parent], make)
        try: 
            os.makedirs(new_dir)
            self.dirs[make] = new_dir
        except: 
            self.add_dir(parent, make)
        return new_dir

    # Given a dataframe, a column, and a suffix, trim the suffix from the column entries and 
    # return the altered dataframe. 
    def truncate(self, dataframe, column, suffix): 
        dataframe = dataframe.copy()
        entries_list = dataframe[column].tolist()
        truncated = [self.truncate_gene(entry, suffix) for entry in entries_list] 
        dataframe[column] = truncated
        return dataframe  
        
    def truncate_gene(self, gene, suffix): 
        try: 
            index = gene.find(suffix)
            return(gene[:index] if index != -1 else gene)
        except: 
            print('Issue detected in trucating gene.') 
            return(gene) 


# In[ ]:





# In[10]:


class NCBI_Synteny(Synteny): 
    def __init__(self, accessions, 
                 run_name, 
                 algs = 'Bfl', 
                 root_directory = './', 
                 proteome_ext = '.faa', 
                 features_ext = '.gff', 
                 genome_ext = '.fna'): 
        self.run_name = run_name
        self.accessions = accessions
        self.dirs = {'root' : root_directory} 
        self.proteome_ext = proteome_ext
        self.features_ext = features_ext
        self.genome_ext = genome_ext 
        self.metadata = self.accessions_table() 
        self.get_algs = self.check_alg_input(algs) 
        self.species_data = {species:dict() for species in self.metadata.index}
        self.build() 
        try: 
            self.make_dir('root', 'ncbi_downloads') 
        except: 
            self.add_dir('root', 'ncbi_downloads')
        print(self.dirs['ncbi_downloads'])
        return 

    def check_alg_input(self, algs): 
        if (algs in self.metadata.index) or (algs == 'Bfl'): 
            return algs 
        elif algs in self.metadata['Accession'].values:
            return self.metadata[self.metadata['Accession']==algs].index[0]
        else: 
            print(algs)
            print(self.metadata['Accession'])
            raise Exception("Please specify an accession or species code for tracing lineage groups, "\
                            "or input 'Bfl' to run with the default bilaterian alg dataset.") 
            
        
    # Make a table containing the metadata from the NCBI assemblies
    def accessions_table(self): 
        datastring = self.get_accession_info()
        metadata = pd.read_csv(StringIO(datastring), sep = '\t') 
        metadata.columns = ['Accession', 'Organism_Name', 'Chromosomes', ]
        binomial_names = metadata['Organism_Name'].tolist()
        species_codes = self.make_species_codes(binomial_names) 
        metadata.index = species_codes
        return metadata

    # With subprocess, pull accession metadata from NCBI using datasets and dataformat
    def get_accession_info(self): 
        try: 
            accessions_string = ','.join(self.accessions)
            datasets_command = ['datasets', 'summary', 'genome', 
                                'accession', accessions_string, '--as-json-lines',]
            # print(' '.join(datasets_command))
            datasets_output = subprocess.run(datasets_command, 
                                             stdout = subprocess.PIPE, text = True, check = True)
            if datasets_output.returncode != 0:
                raise subprocess.CalledProcessError(datasets_output.returncode,
                                                    datasets_command, )
        except Exception as e: 
            raise Exception('Error running datasets in subprocess:') 
            return str(e)
    
        try: 
            dataformat_command = ['dataformat', 'tsv', 'genome', '--fields',
                                  'accession,organism-name,assmstats-total-number-of-chromosomes']
            dataformat_output = subprocess.run(dataformat_command, 
                                               input = datasets_output.stdout, 
                                               stdout = subprocess.PIPE, text = True, check = True)
        except Exception as e: 
            raise Exception('Error formatting datasets in subprocess:') 
            return str(e)
       ## dataframe = pd.DataFrame(dataformat_output.stdout)
        return dataformat_output.stdout

    # Create species codes from the binomial names. Genus species becomes Gsp. 
    # If there are any duplicates, use Gspe instead. 
    def make_species_codes(self, binomial_names): 
        names = [name.split()[0][:1] + name.split()[1][:2] for name in binomial_names]
        if len(names) != len(set(names)): 
            names = [name.split()[0][:2] + name.split()[1][:3] for name in binomial_names]
        return(names) 

    
    ###########
    ### Getting datasets from NCBI 
    ###########
    
    def ncbi_get_datasets(self): 
        with concurrent.futures.ThreadPoolExecutor() as executor: 
            futures = {executor.submit(self.download_and_parse, species) : species
                       for species in self.metadata.index}
            for future in concurrent.futures.as_completed(futures):
                species = futures[future]
                try: 
                    future.result()
                except Exception as e: 
                    print(f"Error downloading or parsing data for {self.metadata.loc[species]['Organism_Name']}.")
                    return e
        if self.get_algs == 'Bfl': 
            self.add_bfl_algs() 
                    
    def download_and_parse(self, species): 
        species_folder = self.make_dir('ncbi_downloads', species) 
        accession = self.metadata.loc[species]['Accession'] 
        datafile =  f"{species_folder}/ncbi_dataset/data/{accession}"
        if os.path.exists(datafile): 
            print(f"Data found for accession {accession}: Skipping {self.metadata.loc[species]['Organism_Name']}.") 
        else: 
            self.ncbi_download(species, species_folder) 
            self.parse_ncbi_zip(species, species_folder) 
        self.arrange_files(species, datafile)
    
    # Downloads datasets from NCBI, but assumes NCBI 'datasets' already installed
    def ncbi_download(self, species, species_folder): 
        accession = self.metadata.loc[species]['Accession']
        filepath = os.path.join(species_folder, 'ncbi_dataset.zip') 
        command = ['datasets',  'download', 
                   'genome', 'accession', accession, 
                   '--include', 'gff3,genome,protein',
                   '--filename', filepath, 
                   '--no-progressbar']
        print(f"Downloading {self.metadata.loc[species]['Organism_Name']} data from accession {accession}. \nThis may take a while...") 
        # print(' '.join(command))
        try: 
            result = subprocess.run(command)
            if result.returncode == 0: 
                print(f"Download complete. Genome {accession} saved to {filepath}.")
            return 
        except: 
            raise Exception(f"Issue encountered while attempting to download accession {accession}.")
    
    def parse_ncbi_zip(self, species, species_folder): 
        filepath = os.path.join(species_folder, 'ncbi_dataset.zip') 
        with zipfile.ZipFile(filepath, 'r') as zip_ref: 
            zip_ref.extractall(species_folder)
        datafile = f"{species_folder}/ncbi_dataset/data/{self.metadata.loc[species]['Accession']}"

    def arrange_files(self, species, datafile): 
        contents = os.listdir(datafile)
        for type, extension in [('proteomes', self.proteome_ext), 
                                ('genomes', self.genome_ext), 
                                ('features', self.features_ext)]: 
            files = [filename for filename in contents if filename.endswith(extension)] 
            
            if len(files) != 1: 
                print(f"Incorrect number of {type} files found for {species} [{len(files)}]") 
            else: 
                identifier = f"{species}_{type}" 
                found_file = os.path.join(datafile, files[0]) 
                new_file_name = os.path.join(datafile, f"{species}{extension}") 
                os.rename(found_file, new_file_name) 
                self.dirs[identifier] = new_file_name

    ###########
    ### Incorporate ALGs from either Bfl or by tracing genes for a given species
    ###########
    def incorporate_algs(self): 
        if self.get_algs == 'Bfl': 
            self.parse_algs('Bfl')
            self.trace_algs('Bfl')
        else: 
            self.trace_chromosomes(self.get_algs)
    
    def add_bfl_algs(self): 
        bfl_path = os.path.join('dependencies', 'Bfl') 
        bfl_files = {'proteomes': 'Bfl_ALGs.fasta', 
                     'genomes': 'Bfl.fna', 
                     'features': 'Bfl_gene_rows.gtf'}
        for type, filename in bfl_files.items(): 
            identifier = 'Bfl_' + type
            self.dirs[identifier] = os.path.join(bfl_path, filename) 
        self.metadata.loc['Bfl'] = ['N/A', 'Branchiostoma floridae', 19]
        self.species_data['Bfl'] = dict() 
        return(bfl_path)

    ############
    ### Building the karyotype file
    ############
    def make_karyotype_files(self): 
        self.build_karyotype()
        columns = ['Chromosome', '1', 'Length', 'SPECIES', '12', '25252']
        labels = ['Chr', 'Start', 'End', 'species', 'size', 'color'] 
        self.clean_karyotype(columns, labels)
        self.write_karyotype() 

    ###########
    ### Executing OrthoFinder
    ###########
    # Make the proteome files, run orthofinder, and incorporate the results 
    def parse_proteomes(self): 
        self.read_proteomes()
        if self.get_algs == 'Bfl': 
            self.proteome_id_trim('Bfl', position = 3) 
        self.proteome_add_species()
        self.write_proteomes()
    
    def execute_orthofinder(self, orthofinder_path = 'orthofinder', threads = None): 
        self.run_orthofinder(orthofinder_path = orthofinder_path, threads = threads)
        self.incorporate_orthofinder()
        self.single_copy_orthologues()
        
    ###########
    ### Building the coordinate files
    ###########
    def build_coordinates(self): 
        for species in self.metadata.index: 
            self.gtf_to_dataframe(species, annotation_type = 'protein_id')
        if self.get_algs == 'Bfl': 
            self.gtf_to_dataframe('Bfl', annotation_type = 'gene_id')
        self.gtf_add_species()
        self.merge_gtf()
        self.incorporate_algs()
        columns = ['annotation', 'sequence', 'start', 'end', 'ALG']
        self.clean_coords(columns)
        self.write_coords()
        
            
   


# In[9]:


def run_NCBI_Synteny(accessions, run_name, algs = 'Bfl', 
                     root_directory = './', orthofinder_path = 'orthofinder', threads = None): 
    synteny = NCBI_Synteny(accessions, run_name, 
                           root_directory = root_directory, algs = algs)
    synteny.ncbi_get_datasets()
    orthofinder_output = os.path.join(synteny.dirs['run_files'], 'orthofinder_output') 
    if not os.path.exists(orthofinder_output): 
        synteny.parse_proteomes()
    synteny.execute_orthofinder(orthofinder_path = orthofinder_path, threads = threads)
    synteny.make_karyotype_files()
    synteny.build_coordinates()
    return(synteny)


# In[ ]:


# class NCBI_Synteny(Synteny): 
#     def __init__(self, accessions, 
#                  species_codes, 
#                  run_name, 
#                  root_directory = './', 
#                  proteome_ext = '.faa', 
#                  features_ext = '.gff', 
#                  genome_ext = '.fna'): 
#         self.species_codes = species_codes
#         self.accessions = dict(zip(species_codes, accessions))
#         self.species_data = {species:dict() for species in species_codes}
#         self.dirs = dict() 
#         self.dirs['root'] = root_directory
#         self.run_name = run_name
#         self.proteome_ext = proteome_ext
#         self.features_ext = features_ext
#         self.genome_ext = genome_ext 
#         self.build()
#         try: 
#             self.make_dir(root, 'ncbi_downloads') 
#         except: 
#             self.add_dir(root, 'ncbi_downloads')
#         return 
        
#     ###########
#     ### Getting datasets from NCBI 
#     ###########
    
#     def ncbi_get_datasets(self): 
#         with concurrent.futures.ThreadPoolExecutor() as executor: 
#             futures = {executor.submit(self.download_and_parse, species) : species
#                        for species in self.accessions.keys()}
#             for future in concurrent.futures.as_completed(futures):
#                 species = futures[future]
#                 try: 
#                     future.result()
#                 except Exception as e: 
#                     print(f"Error downloading or parsing data for species {species}.")
                                       
#     def download_and_parse(self, species): 
#         species_folder = self.make_dir('ncbi_downloads', species) 
#         datafile =  f"{species_folder}/ncbi_dataset/data/{self.accessions[species]}"
#         if os.path.exists(datafile): 
#             print(f"Data found for accession {self.accessions[species]}: Skipping {species}.") 
#         else: 
#             self.ncbi_download(species, species_folder) 
#             self.parse_ncbi_zip(species, species_folder) 
        
#         self.arrange_files(species, datafile)
    
#     # Downloads datasets from NCBI, but assumes NCBI 'datasets' already installed
#     def ncbi_download(self, species, species_folder): 
#         accession = self.accessions[species]
#         filepath = os.path.join(species_folder, 'ncbi_dataset.zip') 
#         command = ['datasets',  'download', 
#                    'genome', 'accession', accession, 
#                    '--include', 'gff3,genome,protein',
#                    '--filename', filepath, 
#                    '--no-progressbar']
#         print(f"Downloading {species} data from accession {accession}. \n This may take a while...") 
#         print(' '.join(command))
#         try: 
#             result = subprocess.run(command)
#             print(f"Download complete. Genome {accession} saved to {filepath}.")
#             return 
#         except: 
#             raise Exception(f"Issue encountered while attempting to download accession {accession}.")
    
#     def parse_ncbi_zip(self, species, species_folder): 
#         filepath = os.path.join(species_folder, 'ncbi_dataset.zip') 
#         with zipfile.ZipFile(filepath, 'r') as zip_ref: 
#             zip_ref.extractall(species_folder)
#         datafile = f"{species_folder}/ncbi_dataset/data/{self.accessions[species]}"

#     def arrange_files(self, species, datafile): 
#         contents = os.listdir(datafile)
#         for type, extension in [('proteomes', self.proteome_ext), 
#                                 ('genomes', self.genome_ext), 
#                                 ('features', self.features_ext)]: 
#             files = [filename for filename in contents if filename.endswith(extension)] 
            
#             if len(files) != 1: 
#                 print(f"Incorrect number of {type} files found for {species} [{len(files)}]") 
#             else: 
#                 identifier = f"{species}_{type}" 
#                 found_file = os.path.join(datafile, files[0]) 
#                 new_file_name = os.path.join(datafile, f"{species}{extension}") 
#                 os.rename(found_file, new_file_name) 
#                 self.dirs[identifier] = new_file_name
        


# In[ ]:


# def check(self, proteome_ext): 
#     self.add_dir('root', 'input_data') 
#     for type, extension in [('proteomes', proteome_ext), 
#                             ('genomes', genome_ext), 
#                             ('features', features_ext)]: 
#         self.add_dir('input_data', type)

#         # print(f"Checking files in {self.dirs[type]}...")
#         contents = os.listdir(self.dirs[type])
#         contents = [filename for filename in contents if filename.endswith(extension)]
#         for species in species_codes: 
#             files = [filename for filename in contents if filename.startswith(species)]
#             # print(files[0])
#             if len(files) != 1: 
#                 print(f"Incorrect number of {type} files found for {species} \n{files}")
#             else: 
#                 identifier = f"{species}_{type}"
#                 self.dirs[identifier] = os.path.join(self.dirs[type], files[0])
#     print('File check complete') 
#     return


# In[ ]:





# In[ ]:




