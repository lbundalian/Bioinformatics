from bs4 import BeautifulSoup
import urllib.request
import json
import csv
import pandas as pd
import requests
import time
import getopt
import json
import re
import wget
import os
import sys


taxa: str = "Hsa"
output_file: str = "ppi.csv"
input_file: str = "genes.txt"
gene_list: str = "DNAH3"
access_key: str = "b29ee574e1fbeb14d7ddb7bcef86afc2"

tax_id: dict = {
    "Mmu": "10090",
    "Hsa": "9606"
}

user_agent = 'Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.0.7) Gecko/2009021910 Firefox/3.0.7'
headers={'User-Agent':user_agent,}


def get_options(argv):

    global taxa
    global output_file
    global input_file


    try:
        opts, args = getopt.getopt(argv,"hi:o:t:",["input=","output=","taxa="])
    except getopt.GetoptError:
        sys.exit(2)

    for opt,arg in opts:
        if opt == '-h':
            print("Available taxa")
            print(taxid)
            sys.exit()
        elif opt in ("-i","--input"):
            input_file = arg
        elif opt in ("-o","--output"):
            output_file = arg
        elif opt in ("-t","--taxa"):
            taxa = arg


get_options(sys.argv[1:])

gene_list = pd.read_csv(f"data/{input_file}",header = None)[0].tolist()
gene_list = "|".join(gene_list)

url = (
    f"https://webservice.thebiogrid.org/interactions?"
    f"searchNames=true&"
    f"geneList={gene_list}&"
    f"includeInteractors=true&"
    f"includeInteractorInteractions=false&"
    f"taxId={tax_id[taxa]}&"
    f"accesskey={access_key}"

)

print("FETCHING FROM BIOGRID!!!!")
request = urllib.request.Request(url,None,headers) #The assembled request
response = urllib.request.urlopen(request)
data = response.read()

print("PARSING RESULTS!!!!")
parsed_result = data.decode().split('\n')
df = pd.DataFrame([x.split('\t') for x in parsed_result]).iloc[:, [7,8]]
df.columns = ['INTERACTOR_A', 'INTERACTOR_B']

print("SAVING OUTPUT!!!!")
df.to_csv(f"output/{output_file}")




