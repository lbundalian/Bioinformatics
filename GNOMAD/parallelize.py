import requests
import sys
import concurrent
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import time
import subprocess

import glob
import os


path: str = sys.argv[2]
out: str = sys.argv[3]
database: str = "/work/users/pz192nijo/Database/gencode.vM29.annotation.gff3.gz"
fmt: str = sys.argv[1]

def load_path(file: str) -> list:
    os.chdir(path)
    files: list = glob.glob(f'*.{fmt}')
    return files
    



def exec_process(file: str) -> str:

    result = "FAIL"

    command = f"htseq-count -f bam -s no -r pos -a 30 {path}/{file} {database} > {out}/{file}.txt"
    try:
        os.system(command)
        #subprocess.run(["htseq-count","bam","no","pos","30", f"{path}/{file}", f"{database}",f"-o {out}/{file}.txt"])
        
    except Exception as e:
        print(e)


    return result



file_list: list = load_path(path)

for i in range(0,len(file_list)):
    with ThreadPoolExecutor(max_workers=len(file_list)) as executor:
        future_to_url = { executor.submit(exec_process, file_list[i])  }
        for future in concurrent.futures.as_completed(future_to_url):
            try:
                data = future.result()
 
            except Exception as e:
                print('Looks like something went wrong:', e)
    
               