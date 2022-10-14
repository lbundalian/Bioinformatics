# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 18:30:49 2022

@author: Linnaeus Bundalian


"""


import requests
import pandas as pd

from bs4 import BeautifulSoup
from bioservices import BioGRID

import urllib.request
import json
import csv
import re
import wget
import os
import sys


# taxid Mmu: 10090
# taxid Hsa: 9606 
taxid: dict = {
        "Mmu": "10090",
        "Hsa": "9606"
    }

degs = pd.read_csv("epigenes.txt",header = None)[0].tolist()

def get_interactors(gene: str,org: str) -> str:

    interactors: str = ""
    try:
        bio_query = BioGRID(query = [f"{gene}"],taxId = taxid["Hsa"])
        bio_interactors = bio_query.biogrid.interactors
        interactors = ','.join([(','.join(map(str, i))) for i in bio_interactors])
    except:
        interactors = ""
    
    return(interactors)



interactors_list = []

for gene in degs:
    i = get_interactors(gene, "Mmu")
    interactors_list.append([gene,i])

df = pd.DataFrame(interactors_list,columns= ['GENE', 'INTERACTORS'])
df.to_csv("epitome.csv")

