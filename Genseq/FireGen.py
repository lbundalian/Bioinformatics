
"""

Author : Linnaeus Bundalian
Date : 25112021
Description : Script to generate dummy multiple sequence alignments based from Jukes-Cantor model for
nucleotide substitution

Features :

    OOP, automatic CLI (fire based)

Last Modified : 28112021

Change logs:
    
    25112021 - Created the file
    26112021 - Implemented Jukes-Cantor Model for mutation
    28112021 - Implemented Jukes-Cantor Model for indels

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


# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSAGenerator(abc.ABC):
    
    @abc.abstractclassmethod
    def save_msa( ):
        pass
    
    @abc.abstractclassmethod
    def generate_seed():
        pass
    
    @abc.abstractclassmethod
    def generate_sequence():
        pass
    
    @abc.abstractclassmethod
    def mutate_sequence():
        pass
    
    @abc.abstractclassmethod
    def print_sequence():
        pass

    
    @abc.abstractclassmethod
    def terminate( ):
        pass
    

class MSAGenerator(IMSAGenerator):
    
    alignments = None
    
    seed_sequence = None
    
    
    nucleotides = ['A','T','G','C']
    stop_codons = ['TAG','TAA','TGA']
    start_codons = ['ATG']

    mutation_prob = []

    use_mutate_prob = False
    use_seed = False
    k = 3
    
    seq_num = 0
    seq_length = 0
    use_seed = False
    per_mutate = 0
    per_shift = 0
    per_introns = 0
    per_stop = 0
    len_introns = 0
    output_file = "output"
    mutated_file = "mutated"
    seed_file = "seed"
    mutation_rate = 0
    indel_rate = 0
    msa = None
    all_codons = []

    
    def __init__(self,n_species, length):

        self.exclude = self.stop_codons + self.start_codons
            
        self.all_codons = [''.join(p) for p in itertools.product(self.nucleotides, repeat=3)] 
        self.codons = [i for i in self.all_codons if i not in self.exclude]

        
        #logging.basicConfig(filename='events.log',level=logging.INFO)
        
        self.seq_num = n_species
        
        self.seq_length = length
        
        self.directory = os.getcwd()
        
        
        
    # =============================================================================
    # Method Name : save_msa
    # Parameters :
    #   msa - the multiple sequence alignment object to be trimmed 
    #   file - the output filename
    #   path -  the working directory of the script 
    # Description :
    # A method used to save a specific msa in the processs
    # =============================================================================

    def save_msa(self,file):
        
        msa = self.msa
        path = self.directory
        
        fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
        for aln in msa:
            SeqIO.write(aln, fasta_file, 'fasta')
        fasta_file.close()

        return self        
    

    def generate_sequence(self,length):

        codon_length = (length/3) - 2
        whole_sequence = ''
        partial_sequence = []
        
        for i in range(int(codon_length)):
            partial_sequence.append(random.choice(self.codons))
        
        whole_sequence = self.start_codons[0] + ''.join(partial_sequence) + ''.join(random.sample(self.stop_codons,1)) 
        
        return whole_sequence

    def generate_seed(self):
    
        self.seed_sequence = self.generate_sequence(self.seq_length)
    
        return self
        
    
    def generate_MSA(self, timer, mutation_rate, indel = True):
    
        
        # self.generate_seed(self.seq_length)
        seq = None
        msa = []
        
        
        for i in range(self.seq_num):
            s =  self.mutate_sequence(self.seed_sequence,timer,mutation_rate)
            if indel:
                s = self.mutate_indel_sequence(s,timer,mutation_rate)
                            
            seq = SeqRecord(Seq(s),
                            id = "Taxon_{0}".format(i+1),
                            name = "Taxon_{0}".format(i+1),
                            description = "")
            msa.append(seq)
                
        self.msa = msa
    
        return self

    
    # number of years?

    def mutate_sequence(self,sequence,t =1, mu = 10):
        
        timer =  t
        seq = list(sequence)
        _seq_length = len(seq)
        
        while timer > 0:
            
            length_time = random.random()
            
            mutation_probability = _seq_length * mu * 0.75
            
            decrement_time = -log(length_time) / mutation_probability
            
            sequence_length = _seq_length

             # index of base to mutate
            base_index = randint(0, sequence_length - 1)
            
            candidate_base = random.choice(self.nucleotides)
            
            # get current base
            current_base = seq[base_index]

            while current_base == candidate_base:
                candidate_base = random.choice(self.nucleotides)

            seq[base_index] = candidate_base
            

            timer -= decrement_time
            
        return ''.join(seq)
         
    
    
    
    def rand_exp(self,_lambda): 
        
        rand_num = -log(random.random()) / _lambda
        
        return rand_num

    
    def poisson_distribution(self,_lambda, T =  1):
        
        times = []
        
        t = self.rand_exp(_lambda)
       
        while t < T: #
            times.append(t)
            t += self.rand_exp(_lambda)
        
        return len(times)
    

    # mutation with indels

    def mutate_indel_sequence(self,sequence,t =1, mu = 10):
        
        timer =  t

        seq = sequence
        _seq_length = len(seq)
        
        
        indel_rate = (_seq_length * mu * timer)/10
        
        n_insert = self.poisson_distribution(indel_rate,timer)
        n_delete = self.poisson_distribution(indel_rate,timer)
        
        
        for n in range(0,n_insert):
            insert_codon = random.choice(self.all_codons)
            insert_index = random.randint(0, len(seq) - 1)
            seq = seq[:insert_index] + insert_codon + seq[-((len(seq)-1)-insert_index):]
        
        
        
        for n in range(0,n_delete):
            delete_index = random.randint(0, len(seq) - 1)
            # print("Deletion expected at {0} for sequence of length {1}".format(delete_index,len(seq)))
            tmp = list(seq)
            tmp.pop(delete_index)
            seq = ''.join(tmp)
            # print("Deleted {0}, new length {1}".format(delete_index,len(seq)))
            
            
        return seq
         



    def print_sequence(self):
        
        for aln in self.msa:
        	
            print(aln.seq)
        
    def terminate(self):
        exit(0)

        
if __name__ == '__main__':
 
    fire.Fire(MSAGenerator)
    
