# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 09092021
Description : Script to remove the part of the sequence which are not conserved across 
the defined species

Last Modified : 22102021

Change logs:
    
    22112021 -  Added step for removing the frameshifted columnn



"""

# =============================================================================
# Imports library 
# =============================================================================

import Bio
import os
import sys
import re
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




threshold = sys.argv[1]
gap = sys.argv[2]
concensus = sys.argv[3]
align1 = sys.argv[4]
species = []


# =============================================================================
# Method Name : remove_gap
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   thr - the % threshold considered for removing the columns with gaps; this will
#         determine the maximum allowable portion of gaps 
# Description :
# A method used to removed sites with gaps among the species of a given multiple
#         (range from 0-100)
# sequence alignment object.
# =============================================================================

def remove_gap(aln,state_cols,kmer = 3):
    
    keep_cols = []
    clean_aln = []

    for i in range(0,len(state_cols),kmer):
        k = 0
        for j in range(kmer):
            k+=state_cols[i+j]
        
        if k == 0:
            for ctr in range(i,i + kmer):
                keep_cols.append(ctr)

    for s in aln:
        sequence = '%s' % ''.join([s.seq[i] for i in keep_cols])
       
        seq = SeqRecord(Seq(sequence),
                        id = s.id,
                        name = s.name,
                        description = s.description)
        clean_aln.append(seq)
    
    return MultipleSeqAlignment(clean_aln)    

# =============================================================================
# Method Name : mask_sequence
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   n - the number of maximum neighbors to mask
# Description :
# A method used to mask n consecutive gaps in a sequence
# =============================================================================

def mask_sequence(aln):
    
    _sequence = []

    for a in aln:
        
        nuc = str(a.seq)
        
        _nucleotides = re.finditer('-{10,}', nuc)
        for n in _nucleotides:
            pos = n.start()
            l = n.end()-n.start()
            nuc = nuc[:pos] + 'X'*l + nuc[pos+l:]	 
            
        nuc = nuc.replace('-','N')
        nuc = nuc.replace('X','-')
        _seq = SeqRecord(Seq(nuc),
                        id = a.id,
                        name = a.name,
                        description = a.description)
        
        _sequence.append(_seq)
    
    return(MultipleSeqAlignment(_sequence))
        

# =============================================================================
# Method Name : create_state_matrix
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   thr - the threshold at which certain block will be flag for discard or keep
# Description :
# A method generate a matrix containing flags for discarding or keeping a block
# =============================================================================

def create_state_matrix(aln, thr, c = '-'):

    clean_aln = []
    state_cols = []
    
    thr = thr/100
    
    for i in range(aln.get_alignment_length()):
        
        col_nucs = [sr.seq[i].upper() for sr in aln]
        counter = Counter(col_nucs)
        # if greater than threshold, discard [1] else keep [0]
        if counter[c]/len(col_nucs) >= thr :
            state_cols.append(1)
        else:
            state_cols.append(0)    

    
    return state_cols

# =============================================================================
# Method Name : create_concensus
# Parameters :
#   aln - the multiple sequence alignment object from which concesus of nucleo-
#   tides will be created 
# Description :
# A method generate a matrix containing concensus for each blocks
# =============================================================================

def create_concensus(aln):

    concensus = []
    
    
    for i in range(aln.get_alignment_length()):
        
        col_nucs = [sr.seq[i].upper() for sr in aln]
        counter = Counter(col_nucs)
        concensus.append(dict(counter))
    
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
    return(MultipleSeqAlignment(_sequence))

# =============================================================================
# Method Name : percent_gap
# Parameters :
#   align - the multiple sequence alignment object to be measured
#   gapthr - the % threshold considered for selecting good quality of alignment among 
#          species
#         (range from 0-100)
# Description :
# A method used to identify species with good alignment in reference to the number
# of gaps.
# =============================================================================

def percent_gap(align,gapthr):
    for a in align:
        nuc_count = Counter(a)
        print(a)
        pergap = (nuc_count['-']/len(a)) * 100  
        if pergap < gapthr:
            _seq = SeqRecord(Seq(a.seq),
                        id = a.id,
                        name = a.name,
                        description = a.description)
        
            species.append(_seq)
        else:
            logging.info("Too much gaps for {0}".format(a.name))
    return(MultipleSeqAlignment(species))

# =============================================================================
# Method Name : find_frameshifts
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
# Description :
# A method to find sequences with frameshifts i.e num_gaps%3 != 0
# =============================================================================

def find_frameshifts(aln,ln = 3,min = 0):
    species = []
    for a in aln:
        nuc = str(a.seq)
        nuc = nuc.replace('N','-')
        matches = re.finditer('-{1,}',nuc)
        mod3 = []
        
        for match in matches:
            mod3.append(len(match.group())%ln)

        if (Counter(mod3)[1] + Counter(mod3)[2]) == min:
            _seq = SeqRecord(Seq(nuc),
                        id = a.id,
                        name = a.name,
                        description = a.description)
            species.append(_seq)
        else:
            logging.info("Frameshift detected in {0}".format(a.name))
    return(MultipleSeqAlignment(species))


def find_frameshifts(aln,ln = 3,min = 0):
    species = []
    for a in aln:
        nuc = str(a.seq)
        nuc = nuc.replace('N','-')
        matches = re.finditer('-{1,}',nuc)
        mod3 = []
        
        for match in matches:
            mod3.append(len(match.group())%ln)

        if (Counter(mod3)[1] + Counter(mod3)[2]) == min:
            _seq = SeqRecord(Seq(nuc),
                        id = a.id,
                        name = a.name,
                        description = a.description)
            species.append(_seq)
        else:
            logging.info("Frameshift detected in {0}".format(a.name))
    return(MultipleSeqAlignment(species))


def remove_frameshifts(aln,ln=3):
    pass



# =============================================================================
# Method Name : matrix_transform
# Parameters :
#   aln - transform the matrix into an array of 3mers 
# Description :
# A method generate a matrix containing flags for discarding or keeping a block
# =============================================================================
def matrix_transform(aln, n = 3):
    msa_matrix = []
    for a in aln:
        _str_aln = str(a.seq)
        str_aln = [_str_aln[i:i+n] for i in range(0, len(_str_aln), n)] 
        msa_matrix.append(str_aln)
    return(np.vstack(msa_matrix))


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
                    print("STOPPED")
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

        return(MultipleSeqAlignment(sequences))
        
        

def extend_MSA(MSA, ext_len):
    
    extended_aln = []
    
    for aln in MSA:
        aln.seq = aln.seq[:-3] + ("X"*ext_len) + aln.seq[-3:]
        _seq = SeqRecord(Seq(aln.seq),
                                    id = aln.id,
                                    name = aln.name,
                                    description = aln.description)
                
        
        extended_aln.append(_seq)
    
    return(MultipleSeqAlignment(extended_aln))


def write_fasta(path,file,MSA):
    
    fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
    for a in MSA:
        
        SeqIO.write(a, fasta_file, 'fasta')
        
    fasta_file.close()

if __name__=='__main__':

    directory = os.getcwd()
    
    logging.basicConfig(filename='events.log',level=logging.INFO)

    debug = False

    # read the alignment
    alignments = AlignIO.read("{0}/{1}".format(directory,align1), "fasta")
        

    if not debug: 
        
        m_aln = len(alignments[0].seq)
        
        if m_aln%3 != 0:
            alignments = extend_MSA(alignments,3-(m_aln%3))
        
        # mask the gaps with bp length < 10
        mask_alignments = mask_sequence(alignments)
        
        write_fasta(directory, "masked", mask_alignments)
        
        # assign if the flags for each blocks [Discard or Keep]
        state_matrix = create_state_matrix(mask_alignments, float(threshold))
        
        # remove the gaps and keep the ones without DISCARD flag - 0
        clean_alignments = remove_gap(mask_alignments,state_matrix)
        
        write_fasta(directory, "clean", clean_alignments)
        
        
        # remove species with %gap
        no_gap_sequences = percent_gap(clean_alignments,float(gap))
        
        write_fasta(directory, "no_gaps", no_gap_sequences)
        
        state_matrix_shifts = create_state_matrix(no_gap_sequences, 10, 'N')
        
        noshift_alignments = remove_gap(no_gap_sequences,state_matrix_shifts,1)
        
        write_fasta(directory, "no_shifts", noshift_alignments)
        
        # remove gaps which is not divisible by 3
        #no_frameshift_sequences = find_frameshifts(no_gap_sequences)   
        
        # check for block wise agreement for the alignments
        # for checking
        con = create_concensus(noshift_alignments)
        matched_matrix = match_concensus(con,noshift_alignments, int(concensus))
        
        no_premature_matrix = search_stop_codon(matched_matrix)
        
        write_fasta(directory, "output", no_premature_matrix)
        
        


