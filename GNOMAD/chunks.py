import requests
import sys
import concurrent
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import time


threads = 29
path: str = "genes"
gnomad: str  = "https://gnomad.broadinstitute.org/api/"

def load_file(path: str,file: str) -> list:
    reads: list = []
    filepath: str = '{0}/{1}'.format(path,file)
    df = pd.read_table(filepath,sep='\t', header=None)
    reads = df
    return reads
    



def get_exac_constraint(gene: str) -> str:

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
            
    except Exception as e:
        print(e.response.text)

    if result == None:
        result = { "pLI": "NOT FOUND","mis_z": "NOT FOUND", "lof_z": "NOT FOUND" }    


    return gene,result['pLI'],result['mis_z'],result['lof_z']



gene_list: list = load_file(path,"epi25_leading.txt")
gene_collection = gene_list[0].values.tolist()
# chunk_list: list = []
# start = int(start)
# offset = int(offset)

# for c in range(start,start+offset):
    # dataset = load_file(path,f"chunk_{c}.csv")
    # chunk_list.append(dataset[0].values.tolist())


n = 29 
chunks = [gene_collection[i:i + n] for i in range(0, len(gene_collection), n)]

gnomad_results: list = []

for i in range(0,len(chunks)):
    tmp: list = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_url = { executor.submit(get_exac_constraint, gene) for gene in chunks[i] }
        for future in concurrent.futures.as_completed(future_to_url):
            try:
                data = future.result()
                tmp.append(data)
                gnomad_results.append(data)
            except Exception as e:
                print('Looks like something went wrong:', e)
    
    pd.DataFrame(tmp).to_csv(f"gnomad_chunks_{i}.csv")            
    time.sleep(60)
               
pd.DataFrame(gnomad_results).to_csv("gnomad_results.csv")