# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 09092021
Description : Script to remove the part of the sequence which are not conserved across 
the defined species

Last Modified : 17092021

"""

# =============================================================================
# Imports library 
# =============================================================================

import Bio
import os
import sys


from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.Align import MultipleSeqAlignment
from collections import Counter




threshold = sys.argv[1]
gap = sys.argv[2]
align1 = sys.argv[3]
species = []



# =============================================================================
# Method Name : remove_gap
# Parameters :
#   aln - the multiple sequence alignment object to be trimmed 
#   thr - the % threshold considered for removing the columns with gaps; this will
#         determine the maximum allowable portion of gaps 
#         (range from 0-100)
# Description :
# A method used to removed sites with gaps among the species of a given multiple
# sequence alignment object.
# =============================================================================


def remove_gap(aln, thr):
    
    clean_aln = []
    keep_cols = []
    
    thr = thr/100
    
    for i in range(aln.get_alignment_length()):
        
        kmer_3 = []

        for j in range(3):
            kmer_3[j] = [sr.seq[i+j].upper() for sr in aln]

        i += 2

        ctr = []   

        #counter = Counter(col_nucs)


        #if counter['-']/len(col_nucs) >= thr :
         #   continue
        for j in range(3):
            
            ctr = Counter(kmer_3[j])

            if ctr['-']/len(kmer_3[j]) >= thr :
                continue

            keep_cols.append(i+j)


    for s in aln:
        sequence = '%s\n' % ''.join([s.seq[i] for i in keep_cols])
       
        seq = SeqRecord(Seq(sequence),
                        id = s.id,
                        name = s.name,
                        description = s.description)
        clean_aln.append(seq)
    
    return MultipleSeqAlignment(clean_aln)


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
        pergap = (nuc_count['-']/len(a)) * 100  
        if pergap < gapthr:
            _seq = SeqRecord(Seq(a.seq),
                        id = a.id,
                        name = a.name,
                        description = a.description)
        
            species.append(_seq)
    return(MultipleSeqAlignment(species))



directory = os.getcwd()

if __name__=='__main__':
    alignments = AlignIO.read("{0}/{1}".format(directory,align1), "fasta")
    clean_alignments = remove_gap(alignments, float(threshold))
    _species = percent_gap(clean_alignments,float(gap))
    for _s in _species:
        print(">%s\n%s" % (_s.id,_s.seq))
        
    #AlignIO.write(_species, "output.fasta", "fasta")
    #os.system('echo {0} > output.fasta'.format(clean_alignments))
    #os.system('{0} {1} {2} {3} > output.fasta'.format('samtool faidx',align1,align2," ".join(species)))

