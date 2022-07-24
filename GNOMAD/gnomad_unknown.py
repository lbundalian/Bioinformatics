import requests
from requests.exceptions import HTTPError
import os
import pandas as pd
import re
import time

path: str = "variants"
gnomad: str  = "https://gnomad.broadinstitute.org/api/"




def get_allele_frequency(url: str, varid: str) -> str:

    allele_frequency: str = ''     

    query: str = """{
      variant(variantId: "var_id", dataset: gnomad_r3) {
        variantId,
        genome{ac,an}
          
        }
    }""".replace("var_id",varid)
    print(varid)
    
    try:
        response: requests.models.Response = requests.post(gnomad, data={'query': query}, timeout = None)
        print(response.status_code)
        if(response.status_code==200):
            print(response.json())
            if response.json().get('errors'):
                allele_frequency = 'NOT_FOUND'
            else:
                allele: dict = response.json()['data']['variant']['genome']
                try:
                    allele_frequency = str(round(allele['ac']/allele['an'],4))
                except Exception as e:
                    allele_frequency = "NA"
        time.sleep(10)
    except Exception as e:
        print(e.response.text)
    
    return allele_frequency

def parse_variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    variant_id: str = None
    variant_id = '{0}-{1}-{2}-{3}'.format(chrom.replace('chr', ''),pos,ref,alt)
    return variant_id

def load_file(path: str) -> list:
    reads: list = []
    for file in os.listdir(path):
        filepath: str = '{0}/{1}'.format(path,file)
        df = pd.read_table(filepath,sep='\t', header=None)
        df = df.loc[:, df.columns!=2]
        df.columns = ['CHROM','POS','AF','REF','ALT','GENE','BIOTYPE','IMPACT','HGVSP']
        reads.append({ 'file':file, 'data':df })
    return reads
        




variant_info: list = load_file(path)

for i in range(0,len(variant_info)):
    variant_info[i]['data']['VarId'] = variant_info[i]['data'].apply(lambda row: parse_variant_id(row['CHROM'],
                    row['POS'],row['REF'],row['ALT']),axis=1)
    variant_info[i]['data']['Gnomad'] = variant_info[i]['data'].apply(lambda row: get_allele_frequency(gnomad,row['VarId']),axis=1)
    variant_info[i]['data'].to_csv(variant_info[i]['file'].replace(".txt",".csv"), index = False)















    
#variant_info[0]['Gnomad'] = variant_info[0].apply(lambda row: get_allele_frequency(gnomad,row['VarId']),axis=1)
#response: requests.models.Response = requests.post(gnomad, data={'query': query}, timeout=None)
#allele: dict = response.json()['data']['variant']['genome']
#allele_frequency: float = round(allele['ac']/allele['an'],4)
