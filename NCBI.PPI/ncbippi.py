# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 18:30:49 2022

@author: Linnaeus Bundalian


"""


import requests
import pandas as pd
import time
import getopt

from bs4 import BeautifulSoup
from bioservices import BioGRID

import urllib.request
import json
import csv
import re
import wget
import os
import sys

delay: float = 0 
taxa: str = "Hsa"
output_file: str = "ppi.csv"
input_file: str = "genes.txt"



taxid: dict = {
        "Mmu": "10090",
        "Hsa": "9606"
    }


def get_options(argv):

    global delay
    global taxa
    global output_file
    global input_file
    
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:d:t:",["input=","output=","delay=","taxa="])
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
        elif opt in ("-d","--delay"):
            delay = float(arg)
        elif opt in ("-t","--taxa"):
            taxa = arg
        
# taxid Mmu: 10090
# taxid Hsa: 9606 


def get_interactors(gene: str,org: str) -> str:

    interactors: str = ""
    try:
        bio_query = BioGRID(query = [f"{gene}"],taxId = taxid[org])
        bio_interactors = bio_query.biogrid.interactors
        interactors = ','.join([(','.join(map(str, i))) for i in bio_interactors])
    except:
        interactors = ""
    
    return(interactors)



interactors_list = []
get_options(sys.argv[1:])

degs = pd.read_csv(input_file,header = None)[0].tolist()


if taxa in taxid: 
    for gene in degs:
        i = get_interactors(gene, taxa)
        time.sleep(delay)
        interactors_list.append([gene,i])

    df = pd.DataFrame(interactors_list,columns= ['GENE', 'INTERACTORS'])
    df.to_csv(output_file)
else:
    print("Please select from the available taxa")