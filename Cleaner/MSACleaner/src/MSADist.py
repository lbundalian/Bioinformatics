# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 26110122
Description : A module for computing

Last Modified : 26012022

Change logs:
    
    26012022 - Created the file

"""

# =============================================================================
# Imports library 
# =============================================================================


import logging
import abc
import distance
import numpy as np


from Bio import pairwise2





# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMISDist(abc.ABC):
    
    
    @abc.abstractclassmethod
    def compute_distance():
        pass
    
    @abc.abstractclassmethod
    def pairwise_distance():
        pass
    
    @abc.abstractclassmethod
    def get_outliers():
        pass
    

class MSADist(IMISDist):
    
    alignments = None
    distances = []
    outliers = []
    
    # =============================================================================
    # Method Name : __init__
    # Parameters :
    #   input_file - file name of the fasta file
    # Description :
    # Constructor method for the class MSADist
    # =============================================================================

    
    def __init__(self, alignments):
        
        logging.basicConfig(filename='events.log',level=logging.INFO)
                
        # assign the alignment
        self.alignments = alignments

    
    
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

    def pairwise_distance(self, reference_id, precision = 3):
        
        self.distances = []
        
        reference = [x for x in self.alignments if (x.id == reference_id )]
        
        
        targets = [x for x in self.alignments if (x.id != reference_id )]
        
        for target in targets:
            result = self.compute_distance(reference[0], target)
            self.distances.append({"Id":target.id,
                                   "Distance": round(result,precision)})
        
                   
        return self
    

    # =============================================================================
    # Method Name : get_outliers
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def get_outliers(self, thr):
        
        self.outliers = []
            
        dist_mean = np.mean([x['Distance'] for x in self.distances])  
        dist_std = np.std([x['Distance'] for x in self.distances])
            
        
        temp = []
            
        for dist in self.distances:
            temp.append({"Id":dist['Id'],
                         "Distance": dist['Distance'],
                         "Score": (dist['Distance'] - dist_mean)/dist_std})
        
        self.outliers = [x for x in temp if abs(x['Score']) > thr]
        
                   
        return self
    
    
    
    
    
    
    

