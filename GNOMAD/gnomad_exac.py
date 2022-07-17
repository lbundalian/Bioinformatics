import requests
from requests.exceptions import HTTPError
import os
import pandas as pd
import re
import time

path: str = "genes"
gnomad: str  = "https://gnomad.broadinstitute.org/api/"




def get_exac_constraint(url: str, gene: str, delay: int) -> str:

    print(gene)
    result = None
    
    query: str = """{
      gene(gene_symbol: "input_gene", reference_genome: GRCh37){
      	gencode_symbol,
        gnomad_constraint{
          pLI,
          mis_z,
          lof_z
        }
      }
    }""".replace("input_gene",gene)
    
    #print(query)
    try:
        response: requests.models.Response = requests.post(gnomad, data={'query': query}, timeout = None)
        #print(response.status_code)
        if(response.status_code==200):
            #print(response.json())
            if response.json().get('errors'):
               result = { "pLI": "NOT FOUND","mis_z": "NOT FOUND", "lof_z": "NOT FOUND" }
            else:
                
                constraint : dict = response.json()['data']['gene']['gnomad_constraint']
                result = constraint    
        time.sleep(delay)
    except Exception as e:
        print(e.response.text)

    if result == None:
        result = { "pLI": "NOT FOUND","mis_z": "NOT FOUND", "lof_z": "NOT FOUND" }    


    return result['pLI'],result['mis_z'],result['lof_z']



def load_file(path: str,file: str) -> list:
    reads: list = []
    filepath: str = '{0}/{1}'.format(path,file)
    df = pd.read_table(filepath,sep='\t', header=None)
    df.columns = ['GENES']
    reads = df 
    return reads
    


    

gene_list: list = load_file(path,"allgenes.txt")





#a = get_exac_constraint(gnomad, "CSMD2",10)
#print(f"{a} {b} {c}")
#print(load_file(path,"epi25.txt"))
#print(gene_list)
gene_list[['pLI','MissenseZ','LoFZ']] = gene_list.apply(lambda row: get_exac_constraint(gnomad,row['GENES'],10),axis=1,result_type = "expand")
gene_list.to_csv("result.csv")



#variant_info[0]['Gnomad'] = variant_info[0].apply(lambda row: get_allele_frequency(gnomad,row['VarId']),axis=1)
#response: requests.models.Response = requests.post(gnomad, data={'query': query}, timeout=None)
#allele: dict = response.json()['data']['variant']['genome']
#allele_frequency: float = round(allele['ac']/allele['an'],4)
