# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 25102021
Description : Script to retrieve CDS data from NCBI. This uses http request libraries 
and scrape data from the server

"""

# =============================================================================
# Imports library 
# =============================================================================

from bs4 import BeautifulSoup
from Bio import Entrez, SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
import urllib.request
import json
import csv
import re
import wget
import os
import sys
import zipfile
import gzip
from distutils.util import strtobool


# =============================================================================
# Initialize variables and accept user input 
# =============================================================================


gene = sys.argv[1]
limit = sys.argv[2]
extract_file = strtobool(sys.argv[3])
remove_files = strtobool(sys.argv[4])

gene_id = ''
fasta_summary = []


user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
headers={'User-Agent':user_agent,}


# =============================================================================
# Methods
# =============================================================================


# =============================================================================
# Method Name : parse_page
# Parameters :
#   url - the link of the website to be read or processed 
# Description :
# Parse a webpage and return the data as text.
# =============================================================================

def parse_page(url):
    request = urllib.request.Request(url,None,headers) #The s request
    response = urllib.request.urlopen(request)
    data = response.read()
    return(data)


# =============================================================================
# Method Name : extract_geneId
# Parameters :
#   html_response - the page parsed by the parse_page method; rendered as text 
# Description :
# Acquire/retrieve the gene id of the gene specied in the automated NCBI genome
# search
# =============================================================================

def extract_geneId(html_response):
    gene_info = BeautifulSoup(html_response, "html.parser").find_all("div", {"class": "ncbi-doc-id"})
    pattern_id = re.compile("[0-9]+")
    gene_id  = pattern_id.search(gene_info[0].text).group(0)
    return(gene_id)


# =============================================================================
# Method Name : retrieve_species
# Parameters :
#   gene_Id - the NCBI id associated with the gene being searched
#   html_response - the page parsed by the parse_page method; rendered as text 
# Description :
# Retrieve a set of species with orthologs of the gene specified
# =============================================================================

def retrieve_species(gene_id, html_response):
    ortholog_info = BeautifulSoup(html_response, "html.parser")
    ortholog_raw = ortholog_info.find_all("script")
    pattern_species = r"\"sci_name\":\ \"[a-zA-Z]+\ [a-zA-Z]+\""
    matches_species = re.finditer(pattern_species, ortholog_raw[11].string, re.MULTILINE)
    species = []
    for m in matches_species:
         species.append(m.group(0)[13:-1])
    return(species)


# =============================================================================
# Method Name : extract_gene_cds
# Parameters :
#   fna_gz - the compressed fasra file downloaded from NCBI
# Description :
# Extracts the coding sequence of the gene specified
# =============================================================================

def extract_gene_cds(fna_gz):
    cds = gzip.open(fna_gz, 'rt')
    fasta = StringIO(cds.read())
    seq = None
    fasta_parsed =  SeqIO.parse(fasta, "fasta")
    gene_pattern = re.compile("gene={0}".format(gene))
    for s in fasta_parsed:
        gene_search = gene_pattern.search(s.description)
        if gene_search != None:
             seq = SeqRecord(Seq(s.seq),
                        id = s.id,
                        name = s.name,
                        description = s.description)
    return(seq)

# =============================================================================
# Method Name : download_CDS
# Parameters :
#   species_list - the list identified by NCBI as species with the gene orthologs
# Description :
# Download the fasta files from NCBI (mammals only)
# =============================================================================

def download_CDS(species_list):
    ctr = 0
    assembly_info = pd.read_csv('Assembly_Reference.csv')[['organism_name','ftp_path','download_path','file_name']]
    for s in species_list:
        ftp_path = assembly_info[assembly_info['organism_name'] == s]['ftp_path']
        if len(ftp_path) != 0:
            url_download = assembly_info[assembly_info['organism_name'] == s]['download_path'].tolist()[0]
            file_to_download = assembly_info[assembly_info['organism_name'] == s]['file_name'].tolist()[0]
            try:
                if not check_file("NCBI//"+file_to_download):
                    result = wget.download(url_download, out = "NCBI//")
                    ctr+=1
                else:
                    ctr+=1
                    
                if limit != 'none' and ctr == int(limit):
                    break
            except:
                print('File not found')


def check_file(file_name):
    return(os.path.isfile(file_name))

def write_to_fasta(sequence):
    fasta_file = open('CDS_{0}.fna'.format(gene),'a+')
    fasta_file.write(">{0}\n".format(str(sequence.description)))
    fasta_file.write(str(sequence.seq)+"\n\n")
    fasta_file.close()

def remove_file(file_name):
    os.remove(file_name)

# =============================================================================
# How it works:
# The script accepts 3 parameters:
# 1. The name of the gene (e.g. CALCR, TRPM3, PDE2A)
# 2. The number of species to be considered in fetching (numeric of range 1 to n; 
#    or none to indicate no limits)
# 3. The indicator if it needs to extract specific gene coding sequence only 
#    ( true : generates gene specific fna file)
# =============================================================================

if __name__=='__main__':

    directory = os.getcwd()
    
    
    
    url_gene = "https://www.ncbi.nlm.nih.gov/homologene/?term={0}".format(gene)
    gene_id = extract_geneId(parse_page(url_gene))
    url_ortholog = "https://www.ncbi.nlm.nih.gov/gene/{0}/ortholog/?scope=40674&term={1}".format(gene_id,gene)
    species = retrieve_species(gene_id, parse_page(url_ortholog))
    
    
    download_CDS(species)
    
    if extract_file:
        ctr = 0
        for files in os.listdir(directory+"//NCBI"):
            if files.endswith('gz'):
                s = extract_gene_cds("NCBI//"+files)
                write_to_fasta(s)
                ctr+=1
                if limit != 'none' and ctr == int(limit):
                        break
            else:
                continue
    
    if remove_files:
        print(remove_files)
        for files in os.listdir(directory+"//NCBI"):
            remove_file(directory+"//NCBI//" + files)
