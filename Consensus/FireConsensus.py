"""

Author : Linnaeus Bundalian
Date : 16122021
Description : Calculate the consensus score based between 2 MSA

Features :

Use case:
    
    python FireGen.py --n_species 1000 --length 10 generate_seed generate_MSA --timer 1 --mutation_rate 0.2 
    --indel True print_sequence

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
import itertools
import random

from math import log
from random import randint

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from itertools import chain
from Bio.Align import MultipleSeqAlignment
from collections import Counter

from numpy.core.fromnumeric import shape
from numpy.core.numeric import normalize_axis_tuple
from numpy.lib.function_base import average

from distutils.util import strtobool

import distance

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSAConsensus(abc.ABC):
    
        
    
    @abc.abstractclassmethod
    def score_consensus():
        pass

    
    @abc.abstractclassmethod
    def terminate( ):
        pass
    
    
class MSAConsensus(IMSAConsensus):

    seed_alignment = None
    new_alignment = None
    directory = None
    input_file = None
    score = 0.0
    
    def __init__(self, seed, fasta):
        
        self.directory = os.getcwd()
        
        self.new_alignment = AlignIO.read("{0}/{1}".format(self.directory,fasta), "fasta")
        self.seed_alignment = AlignIO.read("{0}/{1}".format(self.directory,seed), "fasta")
            
        
    def score_consensus(self):
        
        msa = []
        pairwise_alignment = None
        dummy = []
        
        seed_consensus =  SummaryInfo(self.seed_alignment).dumb_consensus(threshold=0.7)
        target_consensus =  SummaryInfo(self.new_alignment).dumb_consensus(threshold=0.7)
        
        if(len(target_consensus) != len(seed_consensus)):
            pairwise_alignment = pairwise2.align.globalxx(seed_consensus, target_consensus)
            seed_consensus = pairwise_alignment[0].seqA
            target_consensus = pairwise_alignment[0].seqB
        
            
        score = distance.hamming(seed_consensus, target_consensus, True)
        print(score)
        
        
        
    def terminate(self):
        exit(0)

        





if __name__ == '__main__':
 
    fire.Fire(MSAConsensus)
    

    
    