#!/bin/python
# This script will      
#
#
#
# Abin Abraham
# created on: 2018-02-02 22:25:30

import pandas as pd 


INTERSECT_BED_FILE =
MAPPING_DICT_FILE = 
CONSENSUS_FASTA_FILE =  
TE_BED_FILE= "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"


# -----------
# functions
# -----------

def get_intersectbed_file(element, annotation, intersectBedFile_PATH ):
    with open(intersectBedFile_PATH,'r') as f:
        intersectOutput = f.read()
        
    return intersectOutput

def get_conensus_seq_length(CONSENSUS_FASTA_FILE):
    # consensusSeqfile = consensusFastaPath+"consensus_"+element+".fa"
        consensusSeq = SeqIO.read(CONSENSUS_FASTA_FILE,'fasta')
    
    return len(consensusSeq.seq) 

def get_dict_mapping(MAPPING_DICT_FILE):
    # dictname = "all"+element+"_mappedDict.pi"
    allElemMapToCons_dict = pickle.load(open(MAPPING_DICT_FILE, "rb" ))
    
    return allElemMapToCons_dict

# -----------
# main
# ----------- 

#=========== load data into data frame ===========
with open(INTERSECT_BED_FILE,'r') as f:
    intersect_data = f.read()

intersect_data_split =[x.split('\t') for x in intersect_data.split('\n')]
col_name = ["intersect_chr", "intersect_start", "intersect_end", 
            "annotation_score", "TE_chr","TE_start",
            "TE_end", "TE_model","model_start",
            "model_end", "model_length"] 

df = pd.DataFrame(data=out_data[:-1], columns=col_name)
print("Dataframe is loaded.")

#=========== norm index of intesection to genomic instance of element ===========
df["TE_genomicHit_name"] = df["TE_chr"] + ":" + df["TE_start"].map(str) + "-" + df["TE_end"].ma(str)

df[["intersect_start", "intersect_end", "annotation_score", "TE_genomic_start", "TE_genomic_end", "model_start","model_end", "model_length"]]= df[["intersect_start", "intersect_end", "annotation_score", "TE_genomic_start","TE_genomic_end", "model_start", "model_end", "model_length"]].apply(pd.to_numeric, errors='ignore')

df["normToGenomicTE_start"] = df.intersect_start - df.TE_genomic_start 
df["normToGenomicTE_end"] = df.intersect_end - df.TE_genomic_start 
uniq_genomicHits = df.TE_genomicHit_name.unique()
print("Completed mapping annotation coordinates relative to TE start coordinate.")

# =============  Organize Data for Output =============
consensusLength = get_conensus_seq_length(CONSENSUS_FASTA_FILE) # might have to link this to the profile HMM concensus
consArray = np.zeros(consensusLength)
tmp_consArray = np.zeros(consensusLength)
mask_consArray = np.zeros(consensusLength)
print("Done with organizing data.")

# =============  Set Up Mapping Dictionary  =============
# reformat dictionary key 
allMapDict = get_dict_mapping(MAPPING_DICT_FILE)
for k in allMapDict.keys():
    newk = k.split(":")[2] +":"+ k.split(":")[3][:-3]
    renamedMapDict[newk] = allMapDict[k]
print("Mapping Dictionary Loaded")

# =============  Annotate Consensus For Each TE instance =============
annotatedDict = dict()
    for ind,thisElem in enumerate(uniq_genomicHits):
        #filter for one TE
        oneElem_mapped = df.loc[df['TE_genomicHit_name'] == thisElem,["normToGenomicTE_start", "normToGenomicTE_end","annotation_score"]]
        print("currently on: {}, total n = {}".format(str(ind+1), str(len(uniq_genomicHits))))
        
        #get mapping to consensus 
        if renamedMapDict.get(thisElem):
            thisMapDict = renamedMapDict.get(thisElem)
            
            #for this TE, loop through each annotation value (one row per dataframe)
            for ind in range(len(oneElem_mapped)): 
                s,e,v = oneElem_mapped.iloc[ind]
                indexarray = np.arange(s,e)
                #get the indices that match to conensus seq 
                thisconsarrayind = [int(thisMapDict[x]) for x in indexarray if (thisMapDict.get((x)) and int(thisMapDict[x]) >=0)]
                #add annotation to those indices 
                tmp_consArray[thisconsarrayind] = (tmp_consArray[thisconsarrayind] + v)
                #update mask to keep count 
                mask_consArray[thisconsarrayind] += 1

            #divide by mask to get average
            tmp_consArray[mask_consArray>0] = tmp_consArray[mask_consArray>0]/mask_consArray[mask_consArray>0]
            #those indicies without any annotations are converted to nan
            tmp_consArray[mask_consArray==0] = np.nan          
            annotatedDict[thisElem] = tmp_consArray
            
            #reset temporary arrys 
            tmp_consArray = np.zeros(consensusLength)
            mask_consArray = np.zeros(consensusLength)
        else:
            tmp_consArray[mask_consArray==0] = np.nan
            annotatedDict[thisElem] = tmp_consArray
    end = time.time()
    print("Finished annotating consensus sequence in {} minutes.\n".format(round(end-start/60),3))
