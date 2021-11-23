# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 23112021
Description : Script to remove the part of the sequence which are not conserved across 
the defined species

Last Modified : 22112021


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

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.Align import MultipleSeqAlignment
from collections import Counter

from numpy.core.fromnumeric import shape
from numpy.core.numeric import normalize_axis_tuple
from numpy.lib.function_base import average



# =============================================================================
# Initialize variables  
# =============================================================================

colwise_gap_thr = 0
seqwise_gap_thr = 0
concensus_thr = 0

species = []
fasta_file = ''
input_file = ''
output_file = ''


# =============================================================================
# Method Name : get_options
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   thr - the % threshold considered for removing the columns with gaps; this will
#         determine the maximum allowable portion of gaps 
# Description :
# A method used to removed sites with gaps among the species of a given multiple
#         (range from 0-100)
# sequence alignment object.
# =============================================================================

def get_options(argv):

    global colwise_gap_thr
    global seqwise_gap_thr
    global concensus_thr
    global input_file
    global ouput_file
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:b:s:c:",["input=","output=","blockwise=","seqwise=","concensus="])
    except getopt.GetoptError:
        sys.exit(2)
    
    for opt,arg in opts:        
        if opt == '-h':
            print("get lost hahahah") 
            sys.exit()
        elif opt in ("-i","--input"):
            input_file = arg
        elif opt in ("-o","--output"):
            output_file = arg
        elif opt in ("-b","--blockwise"):
            colwise_gap_thr = float(arg)
        elif opt in ("-s", "--seqwise"):
            seqwise_gap_thr = float(arg)
        elif opt in ("-c", "--concensus"):
            concensus_thr = float(arg)


# =============================================================================
# Method Name : save_msa
# Parameters :
#   msa - the multiple sequence alignment object to be trimmed 
#   file - the output filename
#   path -  the working directory of the script 
# Description :
# A method used to save a specific msa in the processs
# =============================================================================

def save_msa(path,file,msa):
    
    fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
    for aln in msa:
        SeqIO.write(aln, fasta_file, 'fasta')
    fasta_file.close()


# =============================================================================
# Method Name : extend_msa
# Parameters :
#   msa - the multiple sequence alignment object to be trimmed 
#   ext_len - the length of the sequence to be inserted 
# Description :
# A method used to extend the MSA so it can be completely divisible by k 
# =============================================================================

def extend_msa(msa, ext_len):
    
    extended_aln = []
    
    for aln in msa:
        aln.seq = aln.seq[:-3] + ("X"*ext_len) + aln.seq[-3:]
        _seq = SeqRecord(Seq(aln.seq),
                                    id = aln.id,
                                    name = aln.name,
                                    description = aln.description)
        extended_aln.append(_seq)
    
    return(MultipleSeqAlignment(extended_aln))



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

def remove_block_gaps(msa,state_matrix,kmer = 3):
        
    keep_cols = []
    no_gaps = []

    # loops through the state matrix in k groups and match to the MSA
    for i in range(0,len(state_matrix),kmer):
        k = 0
        for j in range(kmer):
            k+=state_matrix[i+j]
        
        if k == 0:
            for ctr in range(i,i + kmer):
                keep_cols.append(ctr)

    for alignment in msa:
        sequence = '%s' % ''.join([alignment.seq[i] for i in keep_cols])
        seq = SeqRecord(Seq(sequence),
                        id = alignment.id,
                        name = alignment.name,
                        description = alignment.description)
        no_gaps.append(seq)
    
    return MultipleSeqAlignment(no_gaps)    

# =============================================================================
# Method Name : mask_sequence
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   n - the number of maximum neighbors to mask
# Description :
# A method used to mask n consecutive gaps in a sequence
# =============================================================================

def mask_sequence(msa):
    
    masked_msa = []

    for aln in msa:
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
    
    return MultipleSeqAlignment(masked_msa)
        

# =============================================================================
# Method Name : create_state_matrix
# Parameters :
#   msa - the multiple sequence alignment object to be trimmed 
#   thr - the threshold at which certain block will be flag for discard or keep
#   pattern - the character to whose count will be observed for flagging the blocks
# Description :
# A method generate a matrix containing flags for discarding or keeping a block
# =============================================================================

def create_state_matrix(msa, thr, pattern = '-'):

    state_matrix = []
    
    thr = thr/100
    
    for col in range(msa.get_alignment_length()):
        nuc_cols = [aln.seq[col].upper() for aln in msa]
        counter = Counter(nuc_cols)
        # if greater than threshold, discard [1] else keep [0]
        if counter[pattern]/len(nuc_cols) >= thr :
            state_matrix.append(1)
        else:
            state_matrix.append(0)    

    return state_matrix


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

def remove_seq_gaps(msa, thr):
    
    no_gaps = []
    
    for aln in msa:
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
    
    return MultipleSeqAlignment(no_gaps)


# =============================================================================
# Method Name : create_concensus
# Parameters :
#   aln - the multiple sequence alignment object from which concesus of nucleo-
#   tides will be created 
# Description :
# A method generate a matrix containing concensus for each blocks
# =============================================================================

def create_concensus(msa):

    concensus = []
    
    
    for col in range(msa.get_alignment_length()):
        
        nuc_cols = [aln.seq[col].upper() for aln in msa]
        nuc_counts = Counter(nuc_cols)
        concensus.append(dict(nuc_counts))
    
    return concensus


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

def match_concensus(concensus,aln, thr = 10):

    _sequence = []
    seq_concensus = []
    thr = thr/100

    for a in aln:
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
    
    return MultipleSeqAlignment(_sequence)

# =============================================================================
# Method Name : search_stop_codon
# Parameters :
#   aln - alignment to be observed
# Description :
# A method generate to search for premature stop codon
# =============================================================================
def search_stop_codon(aln, n = 3):
            
        sequences = []

        stop_codons = ['TGA','TAA','TAG']

        for _alignment in aln:
            
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

        return MultipleSeqAlignment(sequences)



# =============================================================================
# Main program 
# =============================================================================

if __name__=='__main__':

    directory = os.getcwd()

    get_options(sys.argv[1:])
    
    logging.basicConfig(filename='events.log',level=logging.INFO)


    debug = False

    # read the alignment
    alignments = AlignIO.read("{0}/{1}".format(directory,input_file), "fasta")
        

    if not debug: 
        
        m_aln = len(alignments[0].seq)
        print(m_aln)
        
        # mask the gaps with bp length < 10
        mask_alignments = mask_sequence(alignments)
        save_msa(directory, "mask", mask_alignments)
 
        
        # remove the masked nucleotides
        state_matrix_shifts = create_state_matrix(mask_alignments, 10, 'N')
        noshift_alignments = remove_block_gaps(mask_alignments,state_matrix_shifts,1)
        save_msa(directory, "noshift", noshift_alignments)
 
    
        m_aln = len(noshift_alignments[0].seq)
        print(m_aln)
        
        
        
        if m_aln%3 != 0:
            noshift_alignments = extend_msa(noshift_alignments,3-(m_aln%3))

            
    
        # assign if the flags for each blocks [Discard or Keep]
        # remove the gaps and keep the ones without DISCARD flag - 0
        state_matrix_gaps = create_state_matrix(noshift_alignments, float(colwise_gap_thr))
        nogap_alignments = remove_block_gaps(noshift_alignments,state_matrix_gaps,3)
        save_msa(directory, "nogaps", nogap_alignments)

        m_aln = len(nogap_alignments[0].seq)
        print(m_aln)
  