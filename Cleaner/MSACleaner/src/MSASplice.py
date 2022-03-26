# =============================================================================
# Imports library 
# =============================================================================

import logging
import abc
import pandas as pd



# =============================================================================
# Creates a formal interface to implement methods on the formal class  
# =============================================================================


class IMSASplice(abc.ABC):
    

    @abc.abstractclassmethod
    def build_splice_table():
        pass
    
    @abc.abstractclassmethod
    def save_splice_table():
        pass

    
    
class MSASplice(IMSASplice):

    # =============================================================================
    # Method Name : __init__
    # Parameters :
    #   input_file - file name of the fasta file
    # Description :
    # Constructor method for the class MSACleaner
    # =============================================================================


    alignments = None
    
    def __init__(self, alignments):
        
        
        
        logging.basicConfig(filename='events.log',level=logging.INFO)
                
        self.alignments = alignments
    
   

    # =============================================================================
    # Method Name : build_splice_stable
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def build_splice_table(self, filename):
        
        
        fo = open(filename, "r")
        
        line = fo.readlines()
        
        results = []
        splice_info = []
        
        index = 0
        ctr = 0
        
        for l in line:
            
            if (l == "== RESULTS ==\n"):
                index = ctr + 1
                
            
            ctr+=1
            results.append(l)
        
        
        for s in range(index,len(results)):
            splice_info.append(results[s].split())
         
            
        splice_df = pd.DataFrame.from_records(splice_info)
        splice_df.columns = ["Gene","Pos","Joint","Sequence","Position","Probability"]
        splice_df = splice_df.drop(['Pos'], axis=1)
        
        self.splice_lookup = splice_df
        
        return self

    # =============================================================================
    # Method Name : save_splice_stable
    # Parameters :
    #   NONE
    # Description :
    # A method used to save the parsed splice table
    # =============================================================================

    def save_splice_table(self):

        self.splice_lookup.to_csv("splice_table.csv") 
        
        return self
        


    # =============================================================================
    # Method Name : check_splice
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================
    def check_splice(self):

        
        
        pass
    
    # =============================================================================
    # Method Name : terminate
    # Parameters :
    #   NONE
    # Description :
    # A method used to terminate the pipeline
    # =============================================================================

    def terminate(self):
        
        exit(0)


    


