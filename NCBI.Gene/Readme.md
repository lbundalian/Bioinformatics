Name : getty.py
Platform : Python v3.6

Directory and file structure:

Make sure that the following file or directories are available:

How the script works:

The script accepts 4 parameters:

Output :

if 3rd parameter is 'True', a fasta file named CDS_{gene_of_interest}.fna is expected

Use case :

python 25102021_NCBI_CDS_Scraping.py CALCR 1 True True

Fetches 1 mammal with ortholog of CALCR
Stores the result to CDS_CALCR.fna
Deletes master fasta file from NCBI directory 