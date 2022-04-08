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
import matplotlib.pyplot as plt

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor




# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSAPhylo(abc.ABC):
    
    
    @abc.abstractclassmethod
    def build_tree():
        pass
    
class MSAPhylo(IMSAPhylo):
    
    alignments = None
    phylo_tree = None
    phylo_lookup = None
    splice_lookup = None
    phylo_parents = {}
    
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
                
        
        return parent_clade
    
    
    def draw_tree(self, mode = 'ascii'):
        
        if (mode == "ascii"):
            Phylo.draw_ascii(self.phylo_tree)        
        elif (mode == "graphics"):
            Phylo.draw(self.phylo_tree)
        else:
            print("Invalid options")


    def save_tree(self):
        fig = plt.figure(figsize=(70, 50), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(self.phylo_tree, axes=axes, do_show=False)     
        plt.savefig('tree.png')


