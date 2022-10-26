# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 25102021
Description : Script to retrieve CDS data from NCBI. This uses http request libraries 
and scrape data from the server

"""

# Imports
from Bio import Entrez
from Bio import SeqIO, AlignIO



# For fetching the mRNA
Entrez.email = "linnaeusbundalian@gmail.com" 
gene = "SIRT7"
taxa = "Homo sapiens"
search_term = f"{gene}[GENE] AND {taxa}[ORGANISM] AND mRNA[Filter] AND RefSeq[Filter]"
print(search_term)
handle = Entrez.read(Entrez.esearch(db = "nucleotide", term = search_term, retmode = "xml"))
genomeIds = handle['IdList']

print(genomeIds)
record = Entrez.efetch(db = "nucleotide", id = genomeIds[0], rettype="gb", retmode="text")
file_out = open("mRNA.gb", "w")
file_out.write(record.read())
