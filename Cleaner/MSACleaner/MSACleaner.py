# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 26110122
Description : An alignment cleaning tool for detecting frameshifts, retained introns and premature stop codons

Last Modified : 26012022

Change logs:
    
    26012022 - Created the file

"""

# =============================================================================
# Imports library 
# =============================================================================


import logging
import abc
import os
import matplotlib.pyplot as plt

from Bio import AlignIO
from Bio import SeqIO
from src.MSADist import MSADist
from src.MSAPhylo import MSAPhylo
from src.MSASplice import MSASplice

# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSACleaner(abc.ABC):

    @abc.abstractclassmethod
    def __init__():
        pass
    
    
    @abc.abstractclassmethod
    def get_outlier_sequence():
        pass
    
    
    @abc.abstractclassmethod
    def save_msa():
        pass
    
    
# Inherits MSADist and MSAPhylo
class MSACleaner(IMSACleaner, MSADist, MSAPhylo, MSASplice):

    directory = None
    input_file = None
    debug = False
    alignments = None
    outlier_sequence = None
    outliers = None
    
    
    
    # =============================================================================
    # Method Name : __init__
    # Parameters :
    #   input_file - file name of the fasta file
    # Description :
    # Constructor method for the class MSACleaner
    # =============================================================================

    
    def __init__(self, input_file, debug):
        
        logging.basicConfig(filename='events.log',level=logging.INFO)
        
        self.directory = os.getcwd()
        
        self.input_file =  input_file
        
        self.debug = debug
        
        # read the alignment
        self.alignments = AlignIO.read("{0}/{1}".format(self.directory,input_file), "fasta")
            
    
    # =============================================================================
    # Method Name : get_outlier_sequence
    # Parameters :
    #   
    # Description :
    # 
    # =============================================================================

    
    def get_outlier_sequence(self):
        
        all_outliers = []
        collect_outliers = []
        
        Species = [taxon.id for taxon in self.alignments]
        
        # perform pairwise hamming distance of nxm alignment; yields n(n-1) 
        # matrix and consolidates the common outliers among the set
        for s in Species:
            
            self.pairwise_distance(s,4)            
            self.get_outliers(3)
            all_outliers.extend(self.outliers)    
        
        
        collect_outliers = list(set([x['Id'] for x in all_outliers]))
  
            
  
        outlier_msa = [s for s in self.alignments if s.id in collect_outliers]        
        
        self.outliers = [x for x in all_outliers if x['Id'] in collect_outliers]
        self.outlier_sequence = outlier_msa
  
        return self
  
    
    def save_msa(self,file):
            
        msa = self.alignments
        path = self.directory    
        
        
        fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
        for aln in msa:
            SeqIO.write(aln, fasta_file, 'fasta')
        
        fasta_file.close()
        
        
     

if __name__ == '__main__':

    cleaner = MSACleaner("aloxe3_cds.fasta", False)
    
    cleaner.build_tree()
    cleaner.save_tree()
    cleaner.get_outlier_sequence()
    cleaner.build_splice_table("splice.out")
    cleaner.save_splice_table()
    print(cleaner.splice_lookup)
    cleaner.save_msa('outlier')
    
    

    
    
    
    
    