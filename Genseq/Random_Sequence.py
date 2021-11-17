"""

Author : Linnaeus Bundalian
Date : 25102021
Description : Script for generating randomized sequence

Update :
    - change mode of accepting console parameters
    - change mode of generating sequences, now uses a seed sequence

"""


# =============================================================================
# Imports library 
# =============================================================================

from Bio import Entrez, SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import unambiguous_dna_by_id

import itertools
import random
import sys
import getopt
import os
import numpy as np
import logging
from distutils.util import strtobool


# =============================================================================
# Initialize variables  
# =============================================================================

# for mutattion the prob to mutate to the other 3 is 80% original 
# is 20% so N1 : 20 N2: 26.67 N3:26.67 N4 : 26.67 



nucleotides = ['A','T','G','C']
stop_codons = ['TAG','TAA','TGA']
start_codons = ['ATG']

mutation_prob = []

use_mutate_prob = False
use_seed = False
k = 3

#seq_length = 3000
#seq_num = 149
#seq_seed = 'ATGATCTAG'


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


seq_seed = ''
exclude = stop_codons + start_codons
codons = [''.join(p) for p in itertools.product(nucleotides, repeat=k)] 
codons = [i for i in codons if i not in exclude]
pseudo_aln = []

def get_options(argv):

    global seq_num
    global seq_length
    global use_seed
    global per_mutate
    global per_shift
    global per_introns
    global len_introns
    global output_file
    global mutated_file
    global per_stop
    
    try:
        opts, args = getopt.getopt(argv,"ho:u:n:m:e:p:s:i:l:a:",["output=","mutated=","n_sequence=","m_length=","use_seed=","p_mutate=","p_shift=","p_introns=","l_introns=","p_stop="])
    except getopt.GetoptError:
        sys.exit(2)
    
    for opt,arg in opts:
        if opt == '-h':
            print("get lost hahahah") 
            sys.exit()
        elif opt in ("-o","--output"):
            output_file = arg
        elif opt in ("-u","--mutated"):
            mutated_file = arg
        elif opt in ("-n", "--n_sequence"):
            seq_num = int(arg)
        elif opt in ("-m", "--m_length"):
            seq_length = int(arg)
        elif opt in ("-e","--use_seed"):
            use_seed = strtobool(arg)
        elif opt in ("-p","--p_mutate"):
            per_mutate = float(arg)
        elif opt in ("-s","--p_shift"):
            per_shift = float(arg)
        elif opt in ("-i","--p_introns"):
            per_introns = float(arg)
        elif opt in ("-l","--l_introns"):
            len_introns = float(arg)
        elif opt in ("-a","--p_stop"):
            per_stop = float(arg)            
            
def generate_sequence(length):
    samples = (length/3) - 2
    randseq = ''
    midseq = []
    for i in range(int(samples)):
        midseq.append(random.choice(codons))
    randseq = start_codons[0] + ''.join(midseq) + ''.join(random.sample(stop_codons,1)) 
    return(randseq)

def generate_seed(length):
    return(generate_sequence(length))

def generate_MSA(num):
    
    seq = None
    seq_seed = generate_seed(seq_length)
    
    if not use_seed:
        for i in range(num):
            s = generate_sequence(seq_length)
            seq = SeqRecord(Seq(s),
                            id = "Species_{0}".format(i+1),
                            name = "Species_{0}".format(i+1),
                            description = "")
            pseudo_aln.append(seq)
    else:
        for i in range(num):
            s = mutate_sequence(seq_seed, per_mutate)
            #print(s)
            seq = SeqRecord(Seq(s),
                            id = "Species_{0}".format(i+1),
                            name = "Species_{0}".format(i+1),
                            description = "")
            pseudo_aln.append(seq)
    return(pseudo_aln)



def write_fasta(path,file,MSA):
    
    fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
    for a in MSA:
        SeqIO.write(a, fasta_file, 'fasta')
        
    fasta_file.close()
       

# for mutating the sequence
# the default probability for mutation
def mutate_sequence(seq, per_mutate):
    
    mutated_seq = ''
    
    n_loc = len(seq)
    for i in range(n_loc)[3:]:
        mutate = random.choices(['Y','N'],[per_mutate, 1-per_mutate])
        if mutate[0] == 'Y':
            tmp = seq[i]
            sequence = list(seq)
            sequence[i] = random.choice([nuc for nuc in nucleotides if nuc != tmp] )
            seq = ''.join(sequence)
    
    mutated_seq = seq

    return(mutated_seq)



# change to actual insertion (addition of sequence) and 
def frameshift(seq, per_shift, p_ins = 0.5):
    
    mutated_seq = ''
    n_loc = int(round(len(seq) * per_shift))
    loc = random.sample([n for n in range(1,len(seq)) if n%3 != 0], n_loc)
    inc_idx = 0
    dec_idx = 0
    
    for i in (range(len(seq))):
        indel_mode = random.choices(['I','D'],[p_ins, 1-p_ins])
        sequence = list(seq)
        if indel_mode[0] == 'D':
            #print("Deleted at {0}".format(i-dec_idx))
            sequence.pop(i-dec_idx)
            seq = ''.join(sequence)
            dec_idx += 1
        else:
            #print("Inserted at {0}".format(i-inc_idx))
            sequence.insert(i+inc_idx,random.choice(nucleotides))
            seq = ''.join(sequence)
            inc_idx += 1
            
    mutated_seq = seq
    return(mutated_seq)



def shift_sequence(aln, per_species):
    
    pseudo_aln = []
    s = ''
    n_species = int(round(len(aln) * per_species))
    idx = random.sample(range(0,len(aln)),n_species)
    for i in range(len(aln)):
        if i in idx:
            logging.info("Frameshift for {0}".format(aln[i].name))
            s = Seq(frameshift(aln[i].seq,0.5,0.5))
        else:
            s = aln[i].seq
        
        pseudo_aln.append(SeqRecord(s,
                            id = "{0}".format(aln[i].id),
                            name = "{0}".format(aln[i].id),
                            description = ""))
    
    return(pseudo_aln)


# this is ok
def premature_stop(seq):
    
    mutated_seq = ''
    kmers = [seq[i:i+3] for i in range(0, len(seq), 3)]
    loc = random.randrange(1,len(kmers)-1,1)
    kmers[loc] = 'TAG'
    mutated_seq = ''.join(str(kmers))
    return(mutated_seq)


def insert_stop(aln, per_stop):
    
    pseudo_aln = []
    s = ''
    n_species = int(round(len(aln) * per_stop))
    idx = random.sample(range(len(aln)),n_species)
    
    for i in range(len(aln)):
        if i in idx:
            logging.info("Premature Stop Codons for {0}".format(aln[i].name))
            s = premature_stop(aln[i].seq)
        else:
            s = aln[i].seq
        
        pseudo_aln.append(SeqRecord(Seq(s),
                            id = "{0}".format(aln[i].id),
                            name = "{0}".format(aln[i].id),
                            description = ""))
    
    return(pseudo_aln)



def generate_introns(length):
    introns = ''.join(np.random.choice(nucleotides, int(length)))
    return(introns)

def retain_introns(aln, per_introns, length_introns):
    
    pseudo_aln = []
    s = ''
    n_species = int(round(len(aln) * per_introns))
    idx = random.sample(range(len(aln)),n_species)
    
    for i in range(len(aln)):
        if i in idx:
            loc = random.choice(range(len(aln[i].seq)))
            logging.info("Intron retention for {0}".format(aln[i].name))
            insert_intron = generate_introns(length_introns)
            s = aln[i].seq[:loc] + insert_intron + aln[i].seq[loc:]
        else:
            s = aln[i].seq
        
        pseudo_aln.append(SeqRecord(s,
                            id = "{0}".format(aln[i].id),
                            name = "{0}".format(aln[i].id),
                            description = ""))
    
    return(pseudo_aln)



if __name__=='__main__':

    directory = os.getcwd()
    get_options(sys.argv[1:])
    logging.basicConfig(filename='{0}_events.log'.format(output_file),level=logging.INFO)
    reference_msa = generate_MSA(seq_num)
    write_fasta(directory,output_file,reference_msa)
    shifted_msa = shift_sequence(reference_msa, per_shift)
    ir_msa = retain_introns(shifted_msa, per_introns, len_introns)
    stop_msa = insert_stop(ir_msa, per_stop)
    write_fasta(directory,mutated_file,stop_msa)