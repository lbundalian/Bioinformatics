Name : NCBI_CDS_Scraping.py
Platform : Python v3.6

Directory and file structure:

Make sure that the following file or directories are available:
1. NCBI
2. Assembly_Reference.csv

How the script works:

The script accepts 4 parameters:
1. The name of the gene (e.g. CALCR, TRPM3, PDE2A)
2. The number of species to be considered in fetching (numeric of range 1 to n; 
   or none to indicate no limits)
3. The indicator if it needs to extract specific gene coding sequence only 
   ( True : generates gene specific fna file)
4. The indicator if user prefer to delete the fasta master file from NCBI 
   ( True : deletes all the files inside NCBI directory)

Output :

if 3rd parameter is 'True', a fasta file named CDS_{gene_of_interest}.fna is expected

Use case :

python 25102021_NCBI_CDS_Scraping.py CALCR 1 True True

Fetches 1 mammal with ortholog of CALCR
Stores the result to CDS_CALCR.fna
Deletes master fasta file from NCBI directory 