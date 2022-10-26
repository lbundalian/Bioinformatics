# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Program: MSc Biomedical Data Science
Date : 22102022
Description : Solution to Task 1 

	a)	Download the sequence of the human SIRT7 gene. Paste the first 20 nucleotides of the gene (complementary strand). 
        Paste the same 20 nucleotides in the main strand. Both sequences should be pasted from 5’ to 3’ (done)
	b)	Paste 100 nucleotides upstream of the gene (in the promoter region)
	c)	Write the positions of the SIRT7 gene in the genomic sequence (the positions should be referenced to the NC_ sequence) 
	d)	Write the positions of the exons in the human SIRT7 gene (done)
	e)	Use the SIRT7 gene sequence and the positions of the exons to make the corresponding SIRT7 transcript (mRNA). Paste the transcript sequence. Indicate the code used to make the transcript (any programming language can be used)
	f)	Write the positions of the 5’UTR, 3’UTR, and CDS of the transcript (done)
	g)	Extract the CDS of the transcript sequence and paste it. Indicate the code used. (done)
	h)	Find the corresponding amino acid sequence from the CDS. Indicate the code used. (done)

"""


# =============================================================================
# Imports library 
# =============================================================================


from Bio import Entrez
from Bio import SeqIO, AlignIO
from Bio import GenBank
from Bio.Seq import Seq

import os
import os.path
import sys
import re
import getopt
import logging
import abc


class IGetty(abc.ABC):

    @abc.abstractclassmethod
    def get_gene_info( ):
        pass
        
    @abc.abstractclassmethod
    def download_genome( ):
        pass

    @abc.abstractclassmethod
    def parse_gene_info( ):
        pass
        
    @abc.abstractclassmethod
    def extract_transcript( ):
        pass
        
    @abc.abstractclassmethod
    def extract_gene( ):
        pass

    @abc.abstractclassmethod
    def print_info( ):
        pass
        
    @abc.abstractclassmethod
    def get_section( ):
        pass

    @abc.abstractclassmethod
    def extract_mrna( ):
        pass
        
    @abc.abstractclassmethod
    def translate_mrna( ):
        pass
        
        
class Getty(IGetty):
    
    # Default
    _email = "linnaeusbundalian@gmail.com" 
    _gene = "SIRT7"
    _taxa = "Homo sapiens"
    _info = {}
    
    _transcript = ""

    def __init__(self,gene:str,taxa:str,email:str):
        
        logging.basicConfig(filename='events.log',level=logging.INFO)
        self._gene = gene
        self._taxa = taxa
        self._email = email
        print("\n")
        print(f"Gene: {self._gene}")
        print(f"Taxa: {self._taxa}")
        print(f"Email: {self._email}")        
        
        self._directory = os.getcwd()
    
    def get_gene_info(self):
        # Getting the gene bank information
        if(os.path.isfile(f"data/{self._gene}.gb")!= True):
            Entrez.email = self._email
            search_term = f"{self._gene}[GENE] AND {self._taxa}[ORGANISM] AND mRNA[Filter] AND RefSeq[Filter]"
            handle = Entrez.read(Entrez.esearch(db = "nucleotide", term = search_term, retmode = "xml"))
            genomeIds = handle['IdList']
            record = Entrez.efetch(db = "nucleotide", id = genomeIds[0], rettype="gb", retmode="text")
            file_out = open(f"data/{self._gene}.gb", "w")
            file_out.write(record.read())
            file_out.close()
            print("\nGene Information Downloaded!!")
        else:
            print(f"\nGENE INFORMATION ALREADY AVAILABLE!")

    def get_genome_info(self):
         
        Entrez.email = self._email
        # Getting gene information
        search_term = f"{self._gene}[GENE] AND {self._taxa}[ORGANISM]"
        handle = Entrez.esearch(db = "gene", term = search_term, retmode = "xml")
        genomeIds = Entrez.read(handle)['IdList']
        record = Entrez.efetch(db = "gene", id = genomeIds[0], rettype="gb", retmode="text")
        info = record.read()
        annot = info.split("\n")[6].split()
        chrom = annot[3]
        start_stop = annot[4].split("..")
        start = start_stop[0][1:]
        stop = start_stop[1][:-1]
        sense = ''.join(char for char in annot[5] if char.isalnum())
        self._info.update({'CHROM': chrom, 'START': start, 'STOP': stop, 'DIRECTION': sense})
 
    def download_genome(self):
        Entrez.email = self._email
        search_term = f"{self._gene}[GENE] AND {self._taxa}[ORGANISM]"
        if(os.path.isfile(f"data/{self._info['CHROM']}.fasta")!= True):
            handle = Entrez.efetch(db = "nucleotide", rettype = "fasta", retmode = "text", id = self._info['CHROM'])
            file_out = open(f"data/{self._info['CHROM']}.fasta", "w")
            print(f"\nDOWNLOADING..........")
            file_out.write(handle.read())
            print(f"\nDOWNLOADED SUCCESSFULLY!!!")
            file_out.close()
        else:
            print(f"\nGENOME FILE IS ALREADY AVAILABLE!")

    def parse_gene_info(self):
       if(os.path.isfile(f"data/{self._gene}.gb")):
            with open(f"data/{self._gene}.gb") as handle:
                for record in GenBank.parse(handle):
                   self._accession = record.accession
                   exons = [x.location for x in record.features if x.key == 'exon']
                   cds = [x.location for x in record.features if x.key == 'CDS']
                   self._info.update({"GENE": self._gene, "ACCESSION": self._accession, "CDS_REGION": cds, "EXONS": exons})
                                      
       else:
            print(f"{self._genbank} NOT FOUND")
        
    def extract_transcript(self):
        
        genome = AlignIO.read("{0}".format(f"data/{self._info['CHROM']}.fasta"), "fasta")
        sequence = str(genome[0].seq)
        sequence = Seq(sequence[int(self._info['START'])-1:int(self._info['STOP'])],)
        if(self._info['DIRECTION']=='complement'):
            self._info.update({'RNA':str(sequence.complement())[::-1]})
        else:
            self._info.update({'RNA':str(sequence)})

        
    def extract_gene(self):
        exon_position = [x.split("..") for x in self._info['EXONS']]
        print(exon_position[0])

    def extract_mrna(self):
        exon_position = [x.split("..") for x in self._info['EXONS']]
        mrna = []
        pre_mrna = self._info['RNA']
        for start,end in exon_position:
            #print(f"{start},{end}")
            mrna.append(pre_mrna[int(start)-1:int(end)])
        mrna = ''.join(mrna)
        self._info.update({'mRNA':mrna})
            
    def extract_cds(self):
        cds_start, cds_stop = self._info['CDS_REGION'][0].split("..")
        cds_start = int(cds_start)
        cds_stop = int(cds_stop)
        cds = self._info['mRNA'][cds_start-1:cds_stop]
        utr5 = self._info['mRNA'][:cds_start-1]
        utr3 = self._info['mRNA'][cds_stop:]
        self._info.update({'CODING': cds,'5UTR': utr5, '3UTR': utr3})
        # print(f"CDS position: {cds_start-1} to {cds_stop}")
        # print(f"5'-UTR position 1 to {cds_start-2}")
        # print(f"3'-UTR position {cds_stop+1} to {len(self._info['mRNA'])}")



    def print_info(self):
        print(f"\nGENE:           {self._info['GENE']}")
        print(f"\nCHROMOSOME:     {self._info['CHROM']}")
        print(f"\nmRNA ACCESSION: {self._info['ACCESSION']}")
        print(f"\nLOCATION:       {self._info['START']},{self._info['STOP']}")
        print(f"\nDIRECTION:      {self._info['DIRECTION']}")
        print(f"\nCDS region:     {self._info['CDS_REGION']}")
        print(f"\nExon location:  {self._info['EXONS']}")
        print(f"\nRNA:            {self._info['RNA']}")
        print(f"\nmRNA:           {self._info['mRNA']}")
        print(f"\nCDS:            {self._info['CODING']}")
        print(f"\n5'UTR:          {self._info['5UTR']}")
        print(f"\n3'UTR:          {self._info['3UTR']}")
        print(f"\nPROTEIN:          {self._info['PROTEIN']}")

    def get_section(self):
        genome = AlignIO.read("{0}".format(f"data/{self._info['CHROM']}.fasta"), "fasta")
        sequence = str(genome[0].seq)
        sequence = Seq(sequence[int(self._info['START'])-1:int(self._info['STOP'])],)
        
        print(f"\nComplementary (20): {str(sequence.complement())[::-1][:20]}")
        print(f"\nMain (20) {str(sequence)[:20]}")
        
        # for the 100 upstream in promoter
        sequence = str(genome[0].seq)
        start = int(self._info['STOP'])+1
        stop = int(self._info['STOP'])+101
        print(f"\nSTART: {start}")
        print(f"STOP: {stop}")
        sequence = Seq(sequence[start:stop],).complement()
        print(f"\nPromoter 100 Upstream: {str(sequence)[::-1]}")
        
    def translate_mrna(self):
        mrna = Seq(self._info['CODING'].replace('T','U'))
        protein = mrna.translate()
        self._info.update({'PROTEIN':protein})
        

        
if __name__ == '__main__':

    getty0 = Getty("SIRT7","Homo sapiens","linnaeusbundalian@gmail")
    getty0.get_genome_info()
    getty0.get_gene_info()
    getty0.parse_gene_info()
    getty0.download_genome()
    getty0.extract_transcript()
    getty0.extract_mrna()
    getty0.extract_cds()
    getty0.translate_mrna()
    getty0.print_info()
    getty0.get_section()
    
    
