# -*- coding: utf-8 -*-
"""

Author : Linnaeus Bundalian
Date : 09092021
Description : Script to remove the part of the sequence which are not conserved across 
the defined species

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



#filename = sys.argv[1]

#parameters = sys.argv

#with open(filename) as file:
#	data = []
#	for line in file:
#		data.append(line[:-1])

#threshold = data[0]
#gap = data[1]
#align = data[2]

threshold = sys.argv[1]
gap = sys.argv[2]
align = sys.argv[3]
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
        
        col_nucs = [sr.seq[i].upper() for sr in aln]
        counter = Counter(col_nucs)

        if counter['-']/len(col_nucs) >= thr :
            continue

        keep_cols.append(i)

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
        if gapthr < pergap:
            species.append(a.id)
    return(species)



directory = os.getcwd()

if __name__=='__main__':
    alignments = AlignIO.read("{0}\\{1}".format(directory,align), "fasta")
    clean_alignments = remove_gap(alignments, float(threshold))
    print(clean_alignments)
    percent_gap(clean_alignments,float(gap))
    print(species)
    os.system('cmd /k {0} {1}'.format('samtool faidx'," ".join(species)))

