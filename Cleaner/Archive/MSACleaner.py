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

import Bio
import os
import sys
import re
import getopt
import numpy as np
import logging
import fire
import abc
import distance
import pandas as pd


from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor




from collections import Counter
from itertools import chain

from numpy.core.fromnumeric import shape
from numpy.core.numeric import normalize_axis_tuple
from numpy.lib.function_base import average
from distutils.util import strtobool


# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSACleaner(abc.ABC):
    
    
    @abc.abstractclassmethod
    def compute_distance():
        pass
    
    @abc.abstractclassmethod
    def pairwise_distance():
        pass
    
    @abc.abstractclassmethod
    def build_tree():
        pass

    @abc.abstractclassmethod
    def build_splice_table():
        pass

    @abc.abstractclassmethod
    def terminate():
        pass
    

class MSACleaner(IMSACleaner):
    
    debug = False
    alignments = None
    directory = ''

    species = []
    fasta_file = ''
    input_file = ''
    output_file = ''
    
    
    distances = []
    phylo_tree = None
    phylo_lookup = None
    splice_lookup = None
    phylo_parents = {}
    
    orf = 0    

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
    # Method Name : compute_distance
    # Parameters :
    #   reference_sequence [Bio.Seq] -  a sequence used to compare with
    #   target_sequence [Bio.Seq] - a sequence under observation
    # Description :
    # Computes the hamming distance between 2 sequences
    # =============================================================================

    def compute_distance(self, reference_sequence, target_sequence):
                
        
        if(len(reference_sequence) != len(target_sequence)):
            pairwise_alignment = pairwise2.align.globalxx(reference_sequence, target_sequence)
            aligned_reference = pairwise_alignment[0].seqA
            aligned_target = pairwise_alignment[0].seqB
        else:
            aligned_reference = reference_sequence
            aligned_target = target_sequence
            
        score = distance.hamming(aligned_reference, aligned_target, True)
        return score     

    # =============================================================================
    # Method Name : pairwise_distance
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def pairwise_distance(self, reference_id):
        
        reference_sequence = None
        
        reference = [x for x in self.alignments if (x.id == reference_id )]
        
        
        targets = [x for x in self.alignments if (x.id != reference_id )]
        
        for target in targets:
            result = self.compute_distance(reference[0], target)
            self.distances.append({"Id":target.id,
                                   "Distance":result})
        
                   
        return self
    
    
    # =============================================================================
    # Method Name : lookup_by_names
    # Parameters :
    #   tree
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    
    def lookup_by_names(self,tree):
        
        names = {}
        
        for clade in tree.find_clades():
            if clade.name:
                if clade.name in names:
                    raise ValueError("Duplicate key: %s" % clade.name)
                names[clade.name] = clade
        return names
    
    
    # =============================================================================
    # Method Name : build_tree
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def build_tree(self):
        
        calculator = DistanceCalculator('identity')
        
        distance_matrix = calculator.get_distance(self.alignments)
        
        # nj or upgm
        constructor = DistanceTreeConstructor(calculator, 'nj')
    
        self.phylo_tree = constructor.build_tree(self.alignments)
        
        self.phylo_lookup = self.lookup_by_names(self.phylo_tree)

        

        for clade in self.phylo_tree.find_clades(order="level"):
            for child in clade:
                self.phylo_parents[child] = clade

            
        return self

    # =============================================================================
    # Method Name : get_parent
    # Parameters :
    #   child [string] - the name of the child node
    # Description :
    # A method returning the parent node of a given child node
    # =============================================================================


    def get_parent(self,child):
        
        child_clade = self.phylo_lookup[child].get_terminals()
        
        parent_clade = self.phylo_parents[child_clade[0]]
                
        print(parent_clade)
        
        return self


    # =============================================================================
    # Method Name : build_splice_stable
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def build_splice_stable(self, filename):
        
        
        fo = open(filename, "r")
        print("Name of the file: ", fo.name)
        
        
        line = fo.readlines()
        
        results = []
        splice_info = []
        
        index = 0
        ctr = 0
        
        for l in line:
            
            if (l == "== RESULTS ==\n"):
                index = ctr + 1
                
            
            ctr+=1
            results.append(l)
        
        
        for s in range(index,len(results)):
            splice_info.results(results[s].split())
         
            
        splice_df = pd.DataFrame.from_records(splice_info)
        splice_df.columns = ["Gene","Pos","Joint","Sequence","Position","Probability"]
        splice_df = splice_df.drop(['Pos'], axis=1)
        
        self.splice_lookup = splice_df
        
        return self


    # =============================================================================
    # Method Name : check_splice
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================
    def check_splice(self):

        
        
        pass
    
    # =============================================================================
    # Method Name : terminate
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def terminate(self):
        
        exit(0)

        
if __name__ == '__main__':
 
    fire.Fire(MSACleaner)
    

    


