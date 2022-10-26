# genbank

from Bio import Entrez
from Bio import SeqIO, AlignIO
from Bio import GenBank

from Bio import GenBank
with open("data/SIRT7.gb") as handle:
    for record in GenBank.parse(handle):
        print(record.accession)
        #print(record.features[3].key)
        exons = [x.location for x in record.features if x.key == 'exon']
        CDS = [x.location for x in record.features if x.key == 'CDS']
        print(exons)
        print(CDS)