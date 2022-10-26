# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 25112021
Description : Script to remove the part of the sequence which are not conserved across 
the defined species

Last Modified : 25112021

Change logs:
    
    25112021 - Created the file

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


class IMSACleaner(abc.ABC):
    
    @abc.abstractclassmethod
    def save_msa( ):
        pass
    
    @abc.abstractclassmethod
    def extend_msa( ):
        pass
    
    @abc.abstractclassmethod
    def find_orf( ):
        pass
    
    @abc.abstractclassmethod
    def remove_block_gaps( ):
        pass
    
    @abc.abstractclassmethod
    def mask_sequence( ):
        pass
    
    @abc.abstractclassmethod
    def create_state_matrix( ):
        pass
    
    @abc.abstractclassmethod
    def remove_seq_gaps( ):
        pass
    
    @abc.abstractclassmethod
    def create_concensus( ):
        pass
    
    @abc.abstractclassmethod
    def match_concensus( ):
        pass

    @abc.abstractclassmethod
    def search_stop_codon( ):
        pass
    
    @abc.abstractclassmethod
    def margin_trim( ):
        pass
    
    @abc.abstractclassmethod
    def find_frameshifts( ):
        pass

    @abc.abstractclassmethod
    def default_pipeline( ):
        pass

    
    @abc.abstractclassmethod
    def terminate( ):
        pass
    

class MSACleaner(IMSACleaner):
    
    colwise_gap_thr = 0
    strict_shift_thr = 95
    strict_gap_thr = 90
    seqwise_gap_thr = 0
    concensus_thr = 0
    debug = False
    alignments = None
    directory = ''

    species = []
    fasta_file = ''
    input_file = ''
    output_file = ''
    orf = 0
    
    state_matrix = []
    
    
    def __init__(self, input_file):
        
        logging.basicConfig(filename='events.log',level=logging.INFO)
        
        self.directory = os.getcwd()
        
        self.input_file =  input_file
        
        # read the alignment
        self.alignments = AlignIO.read("{0}/{1}".format(self.directory,input_file), "fasta")
        
        
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
        
        msa = self.alignments
        path = self.directory
        
        fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
        for aln in msa:
            SeqIO.write(aln, fasta_file, 'fasta')
        fasta_file.close()

        return self        
    
    
    # =============================================================================
    # Method Name : extend_msa
    # Parameters :
    #   msa - the multiple sequence alignment object to be trimmed 
    #   ext_len - the length of the sequence to be inserted 
    # Description :
    # A method used to extend the MSA so it can be completely divisible by k 
    # =============================================================================

    def extend_msa(self, ext_len):
        
        extended_aln = []
        
        for aln in self.alignments:
            aln.seq = aln.seq[:-3] + ("X"*ext_len) + aln.seq[-3:]
            _seq = SeqRecord(Seq(aln.seq),
                                        id = aln.id,
                                        name = aln.name,
                                        description = aln.description)
            extended_aln.append(_seq)
        
        self.alignments = MultipleSeqAlignment(extended_aln)
        
        return self

    
    # =============================================================================
    # Method Name : find_orf
    # Parameters :
    #   aln - the multiple sequence alignment object from to concesus of nucleo-
    # Description :
    # A method remove non conforming sequence
    # =============================================================================

    def find_orf(self):

        positions = []
        position_counts = []
        
        for aln in self.alignments:
            nuc = str(aln.seq)
            _nucleotides = re.finditer('ATG+', nuc)
            tmp = []
            for n in _nucleotides:
                pos = n.start()
                tmp.append(pos)
            positions.append(tmp)
            
        positions_combined = list(chain.from_iterable(positions))
        position_counts = Counter(positions_combined)
        orf,n = position_counts.most_common(1)[0]
        
        
        return orf
    
    
    # =============================================================================
    # Method Name : remove_gap
    # Parameters :
    #   msa - the multiple sequence alignment object to be trimmed 
    #   state_matrix - a matrix containing the discard or keep flags of the blocks
    #   kmer - the number of nucleotides in a search group or pattern   
    # Description :
    # A method used to removed sites with gaps among the species of a given multiple
    #         (range from 0-100)
    # sequence alignment object.
    # =============================================================================
    
    def remove_block_gaps(self,kmer = 3):
            
        keep_cols = []
        no_gaps = []
        
        state_matrix = self.state_matrix
    
        # loops through the state matrix in k groups and match to the MSA
        for i in range(0,len(state_matrix),kmer):
            k = 0
            for j in range(kmer):
                k+=state_matrix[i+j]
            
            if k == 0:
                for ctr in range(i,i + kmer):
                    keep_cols.append(ctr)
    
        for alignment in self.alignments:
            sequence = '%s' % ''.join([alignment.seq[i] for i in keep_cols])
            seq = SeqRecord(Seq(sequence),
                            id = alignment.id,
                            name = alignment.name,
                            description = alignment.description)
            no_gaps.append(seq)
        
        self.alignments = MultipleSeqAlignment(no_gaps)
        
        return self    
    
    
    # =============================================================================
    # Method Name : mask_sequence
    # Parameters :
    #   aln - the multiple sequence alignment object to be trimmed 
    #   n - the number of maximum neighbors to mask
    # Description :
    # A method used to mask n consecutive gaps in a sequence
    # =============================================================================
    
    def mask_sequence(self):
        
        masked_msa = []
    
        for aln in self.alignments:
            nuc = str(aln.seq)
            _nucleotides = re.finditer('-{10,}', nuc)
            for n in _nucleotides:
                pos = n.start()
                l = n.end()-n.start()
                nuc = nuc[:pos] + 'X'*l + nuc[pos+l:]	            
            nuc = nuc.replace('-','N')
            nuc = nuc.replace('X','-')
            _seq = SeqRecord(Seq(nuc),
                            id = aln.id,
                            name = aln.name,
                            description = aln.description)
            masked_msa.append(_seq)
        
        self.alignments = MultipleSeqAlignment(masked_msa)
        
        return self
        
    # =============================================================================
    # Method Name : create_state_matrix
    # Parameters :
    #   msa - the multiple sequence alignment object to be trimmed 
    #   thr - the threshold at which certain block will be flag for discard or keep
    #   pattern - the character to whose count will be observed for flagging the blocks
    # Description :
    # A method generate a matrix containing flags for discarding or keeping a block
    # =============================================================================
    
    def create_state_matrix(self, thr, pattern = '-'):
    
        state_matrix = []
        
        thr = thr/100
        
        for col in range(self.alignments.get_alignment_length()):
            nuc_cols = [aln.seq[col].upper() for aln in self.alignments]
            counter = Counter(nuc_cols)
            # if greater than threshold, discard [1] else keep [0]
            if counter[pattern]/len(nuc_cols) >= thr :
                state_matrix.append(1)
            else:
                state_matrix.append(0)    
    
        self.state_matrix = state_matrix    
    
        return self



    # =============================================================================
    # Method Name : remove_seqgaps
    # Parameters :
    #   align - the multiple sequence alignment object to be measured
    #   gapthr - the % threshold considered for selecting good quality of alignment among 
    #          species
    #         (range from 0-100)
    # Description :
    # A method used to identify species with good alignment in reference to the number
    # of gaps.
    # =============================================================================

    def remove_seq_gaps(self, thr):
        
        no_gaps = []
        
        for aln in self.alignments:
            nuc_counts = Counter(aln)        
            percent_gap = (nuc_counts['-']/len(aln)) * 100
            if percent_gap <= thr:
                _seq = SeqRecord(Seq(aln.seq),
                            id = aln.id,
                            name = aln.name,
                            description = aln.description)
                no_gaps.append(_seq)
            else:
                logging.info("Too much gaps for {0}".format(aln.name))
        
        self.alignments =  MultipleSeqAlignment(no_gaps)
        
        return self

    # =============================================================================
    # Method Name : create_concensus
    # Parameters :
    #   aln - the multiple sequence alignment object from which concesus of nucleo-
    #   tides will be created 
    # Description :
    # A method generate a matrix containing concensus for each blocks
    # =============================================================================

    def create_concensus(self):

        concensus = []
        
        
        for col in range(self.alignments.get_alignment_length()):
            
            nuc_cols = [aln.seq[col].upper() for aln in self.alignments]
            nuc_counts = Counter(nuc_cols)
            concensus.append(dict(nuc_counts))
        
        self.concensus = concensus
        
        return self

    # =============================================================================
    # Method Name : match_concensus
    # Parameters :
    #   aln - the multiple sequence alignment object from to concesus of nucleo-
    #   tides will be matched
    #   concensus -  the amount of each nucleotides in a block in batches
    #   thr -  the allowable % portion of the sequence which can be considered as 
    #   non conforming nucleotide across a block
    # Description :
    # A method remove non conforming sequence
    # =============================================================================

    def match_concensus(self, thr = 10):

        _sequence = []
        seq_concensus = []
        thr = thr/100
        concensus = self.concensus

        for a in self.alignments:
            for _a in range(len(a.seq)):
                
                # and the portion at which the nucleotide constitute should be larger than the average or larger than 20%

                if concensus[_a][a.seq[_a]] > concensus[_a][min(concensus[_a], key = concensus[_a].get)]:
                    # good
                    seq_concensus.append(1)
                else :
                    # bad
                    seq_concensus.append(0)
        
            ctr = dict(Counter(seq_concensus))
            per_ctr = ctr[0]/sum(ctr.values())
            if( per_ctr < thr):
                _seq = SeqRecord(Seq(a.seq),
                            id = a.id,
                            name = a.name,
                            description = a.description)
            
            
                _sequence.append(_seq)
            else:
                logging.info("Too high percent mismatch in concensus in {0} - {1:.4g}%".format(a.name,per_ctr*100))
        
        self.alignments = MultipleSeqAlignment(_sequence)
        
        return self

    # =============================================================================
    # Method Name : search_stop_codon
    # Parameters :
    #   aln - alignment to be observed
    # Description :
    # A method generate to search for premature stop codon
    # =============================================================================
    def search_stop_codon(self, n = 3):
                
            sequences = []

            stop_codons = ['TGA','TAA','TAG']

            for _alignment in self.alignments:
                
                premature = False
                
                for index in range(0, len(_alignment.seq), 3):
                    codon = _alignment.seq[index:index+3]
                    if codon in stop_codons and index < (len(_alignment.seq)/3):
                        logging.info("Premature stop codon in {0} Block : {1}".format(_alignment.name, index))
                        premature = True
                        break

                

                if not premature:
                    
                    _seq = SeqRecord(Seq(_alignment.seq),
                                        id = _alignment.id,
                                        name = _alignment.name,
                                        description = _alignment.description)
                    
            
                    sequences.append(_seq)
          
                else:
                    
                    continue

            self.alignments = MultipleSeqAlignment(sequences)
            
            return self 


    # =============================================================================
    # Method Name : margin_trim
    # Parameters :
    #   msa - alignment to be observed
    #   start - 
    # Description :
    # A method generate to search for premature stop codon
    # =============================================================================
    def margin_trim(self, start = None, end = None):

        trimmed_msa = []    

        for aln in self.alignments:
            _seq = SeqRecord(Seq(aln.seq[start:end]),
                            id = aln.id,
                            name = aln.name,
                            description = aln.description)
            trimmed_msa.append(_seq)

        self.alignments = trimmed_msa
       
        return self


    # =============================================================================
    # Method Name : find_frameshifts
    # Parameters :
    #   aln - the multiple sequence alignment object to be trimmed 
    # Description :
    # A method to find sequences with frameshifts i.e num_gaps%3 != 0
    # =============================================================================

    def find_frameshifts(self):
        
        noshifts = []
        
        
        for aln in self.alignments:
            n_position = []        
            nuc = str(aln.seq)
            matches = re.finditer('N{1,}',nuc)
            for match in matches:
                n_position.append(match.start())
            
            
            
            if len(n_position) == 0:
                _seq = SeqRecord(Seq(nuc),
                            id = aln.id,
                            name = aln.name,
                            description = aln.description)
                noshifts.append(_seq)
            else:
                logging.info("Frameshift detected in {0}".format(aln.name))
        
        self.alignments = MultipleSeqAlignment(noshifts)

        return self

    # =============================================================================
    # Method Name : default_pipeline
    # Parameters :
    #   aln - the multiple sequence alignment object to be trimmed 
    # Description :
    # A method to find sequences with frameshifts i.e num_gaps%3 != 0
    # =============================================================================
    def default_pipeline(self, strict_shift_thr, strict_gap_thr, colwise_gap_thr, seqwise_gap_thr, concensus_thr):
        
        # find the ORF
        self.margin_trim(start = self.find_orf())
        
        # mask the gaps with bp length < 10
        self.mask_sequence()
        
        # remove the masked nucleotides
        self.create_state_matrix(strict_shift_thr,'N')
        self.remove_block_gaps(1)
        
        # High gap percentage, this will more likely remove insertions
        self.create_state_matrix(strict_gap_thr)
        self.remove_block_gaps(1)
       
        
        # assign if the flags for each blocks [Discard or Keep]
        # remove the gaps and keep the ones without DISCARD flag - 0
        self.create_state_matrix(float(colwise_gap_thr))
        self.remove_block_gaps(3)
       
        # remove gappy sequences
        self.remove_seq_gaps(float(seqwise_gap_thr))
        
        
        self.find_frameshifts()
       
        self.create_concensus()
        
        self.match_concensus(concensus_thr)

        self.save_msa("pipeline_output")
        
        self.search_stop_codon()
        self.save_msa("pipeline_output_1")
        # save_msa(directory + "/Logs", "nostop", nopremature_alignments)


    
    def terminate(self):
        exit(0)

    def print_sequence(self):
        for aln in self.alignments:
            print(aln.seq)
        
if __name__ == '__main__':
 
    fire.Fire(MSACleaner)
    
