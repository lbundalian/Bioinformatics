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
from Bio import GenBank


# For fetching the mRNA
# Entrez.email = "linnaeusbundalian@gmail.com" 
# gene = "SIRT7"
# taxa = "Homo sapiens"
# search_term = f"{gene}[GENE] AND {taxa}[ORGANISM]"#  AND mRNA[Filter] AND RefSeq[Filter]"
# print(search_term)
# handle = Entrez.read(Entrez.esearch(db = "nucleotide", term = search_term, retmode = "xml"))
# genomeIds = handle['IdList']

# print(genomeIds)
# record = Entrez.efetch(db = "nucleotide", id = genomeIds[0], rettype="fasta", retmode="text")
# file_out = open("transcript.fasta", "w")
# file_out.write(record.read())
# file_out.close()


#for fetching the whole genome
Entrez.email = "linnaeusbundalian@gmail.com" 
gene = "SIRT7"
taxa = "Homo sapiens"
search_term = f"{gene}[GENE] AND {taxa}[ORGANISM]"#  AND mRNA[Filter] AND RefSeq[Filter]"
print(search_term)
handle = Entrez.esearch(db = "gene", term = search_term, retmode = "xml")

genomeIds = Entrez.read(handle)['IdList']

print(genomeIds)
record = Entrez.efetch(db = "gene", id = genomeIds[0], rettype="gb", retmode="text")
info = record.read()

print(info)

annot = info.split("\n")[6].split()
print(annot)
chrom = annot[3]
start_stop = annot[4].split("..")
start = start_stop[0][1:]
stop = start_stop[1][:-1]
print(f"Chromid: {chrom} Start:{start} Stop: {stop}")


# Entrez.email = "linnaeusbundalian@gmail.com" 
# search_term = f"{gene}[GENE] AND {taxa}[ORGANISM]"#  AND mRNA[Filter] AND RefSeq[Filter]"
# print(search_term)
# handle = Entrez.efetch(db = "nucleotide", rettype = "fasta", retmode = "text", id = chrom)


# #record = Entrez.efetch(db = "nucleotide", id = genomeIds[0], rettype="fasta", retmode="text")

# file_out = open("transcript.fasta", "w")
# print("START")
# file_out.write(handle.read())
# print("END")
# file_out.close()


# Import Libraries
from Bio.Seq import Seq

genome = AlignIO.read("{0}".format("transcript.fasta"), "fasta")
sequence = str(genome[0].seq)
sequence = Seq(sequence[int(81911939)-1:int(81918176)],)
#print(sequence)
print(str(sequence.complement())[::-1])


