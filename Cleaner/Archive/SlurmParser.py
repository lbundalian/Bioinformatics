# Open a file

import pandas as pd

fo = open("slurm.out", "r")
print("Name of the file: ", fo.name)

# Assuming file has following 5 lines
# This is 1st line
# This is 2nd line
# This is 3rd line
# This is 4th line
# This is 5th line

line = fo.readlines()

slurm = []
splice_info = []

index = 0
ctr = 0
for l in line:
    
    if (l == "== RESULTS ==\n"):
        index = ctr + 1
        
    
    ctr+=1
    slurm.append(l)


for s in range(index,len(slurm)):
    splice_info.append(slurm[s].split())
 
    
splice_df = pd.DataFrame.from_records(splice_info)
splice_df.columns = ["Gene","Pos","Joint","Sequence","Position","Probability"]
splice_df = splice_df.drop(['Pos'], axis=1)



fo.close()