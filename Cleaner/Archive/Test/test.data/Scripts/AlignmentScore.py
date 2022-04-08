from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from itertools import chain
from Bio.Align import MultipleSeqAlignment
from collections import Counter



from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62, FileMatrix
from pymsa.util.fasta import print_alignment


import os
import sys
import re
import getopt


fasta_file = ''

def run_all_scores(sequences: list) -> None:
    aligned_sequences = list(pair[1] for pair in sequences)
    sequences_id = list(pair[0] for pair in sequences)

    msa = MSA(aligned_sequences, sequences_id)
    

    # Percentage of non-gaps and totally conserved columns
    non_gaps = PercentageOfNonGaps(msa)
    totally_conserved_columns = PercentageOfTotallyConservedColumns(msa)

    percentage = non_gaps.compute()
    print("Percentage of non-gaps: {0} %".format(percentage))

    conserved = totally_conserved_columns.compute()
    print("Percentage of totally conserved columns: {0}".format(conserved))

    # Entropy
    value = Entropy(msa).compute()
    print("Entropy score: {0}".format(value))

    # Sum of pairs
    value = SumOfPairs(msa, Blosum62()).compute()
    print("Sum of Pairs score (Blosum62): {0}".format(value))

    value = SumOfPairs(msa, PAM250()).compute()
    print("Sum of Pairs score (PAM250): {0}".format(value))

#    value = SumOfPairs(msa, FileMatrix('PAM380.txt')).compute()
 #   print("Sum of Pairs score (PAM380): {0}".format(value))

    # Star
    value = Star(msa, Blosum62()).compute()
    print("Star score (Blosum62): {0}".format(value))

    value = Star(msa, PAM250()).compute()
    print("Star score (PAM250): {0}".format(value))
    
if __name__ == '__main__':
    
    directory = os.getcwd()
    alignments = None

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:",["input="])
    except getopt.GetoptError:
        sys.exit(2)

    for opt,arg in opts:       
        if opt == '-h':
            print("get lost hahahah") 
            sys.exit()
        elif opt in ("-i","--input"):
            fasta_file = arg

    alignments = AlignIO.read("{0}/{1}".format(directory,fasta_file), "fasta")

    sequences = [] 
    
    for aln in alignments:
        sequences.append((aln.id,str(aln.seq)))
    
    print("Length of sequence: {0}".format(len(str(alignments[0].seq))))
    run_all_scores(sequences)
    
