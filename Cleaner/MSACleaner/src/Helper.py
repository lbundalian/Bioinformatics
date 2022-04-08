import sys
from Bio import SeqIO


def progress(count, total, status=''):
    sys.stdout.write('\r')
    
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = u"\u2588" * filled_len + ' ' * (bar_len - filled_len)

    sys.stdout.write('|%s| %s%s' % (bar, percents, '%'))
    sys.stdout.flush() 


def zscore(value, mean, std):
    zscore = 0
    zscore = (value - mean)/std
    return zscore 

def save_msa(alignments,path,file):
        
    msa = alignments
        
    fasta_file = open('{0}//{1}.fasta'.format(path,file),'w+')
    for aln in msa:
        SeqIO.write(aln, fasta_file, 'fasta')
    
    fasta_file.close()

    