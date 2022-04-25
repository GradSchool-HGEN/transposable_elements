#!/bin/python
# This script will plots the number of bases in aligned consensus to number of bases in aligned genomic hit. Then colors each TE based on whether the element overlaps a enhancer.abs   
#
#
#
# Abin Abraham
# created on: ...2017-12-19 13:19:09s

import numpy as np 
import pandas as pd 
import subprocess 
import matplotlib.pyplot as plt 


element = "MER20"
TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
annotate_enhancer_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/all_facet-organ_enhancers_coordinates-formated.tsv"
filterForElement = subprocess.Popen(['grep', '-w', element, TEbedfile_PATH], stdout=subprocess.PIPE, universal_newlines=True)
intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-loj',  '-a', 'stdin','-b', annotate_enhancer_PATH, '-sorted'], stdin = filterForElement.stdout, universal_newlines=True)  
    


out_data =[x.split('\t') for x in intersectOutput.split('\n')]
col_name = ["TEchr", "TEstart", "TEend", "TE", "model_start","model_end","model_length", "intersect_chr","intersect_start", "intersect_end" , "intersect_value"]
df = pd.DataFrame(data=out_data[:-1], columns=col_name)
df = df.apply(pd.to_numeric, errors='ignore')

df['genomicNumBases'] = df.TEend - df.TEstart 
df['alignConsensNumBases'] = df.model_end - df.model_start
df['deviationFromConsensus'] = abs(df.genomicNumBases-df.alignConsensNumBases)
df['per_deviationFromConsensus'] = abs(df.genomicNumBases-df.alignConsensNumBases)/df.model_length

df['enhancerOverlap'] = 0
df.loc[df.intersect_chr != ".", 'enhancerOverlap'] = 1
df.loc[df.intersect_value == 1]

# plot 

#=========== scatter plots ===========
plt.figure(0)
plt.scatter(df.loc[df.enhancerOverlap == 0, 'genomicNumBases'], df.loc[df.enhancerOverlap == 0, 'alignConsensNumBases'], color='b',  alpha=0.2,  marker=".", s =2, label="No Enhancer Overlap (n= "+str(df.loc[df.enhancerOverlap == 0, 'genomicNumBases'].shape[0])+")")
plt.scatter(df.loc[df.enhancerOverlap == 1, 'genomicNumBases'], df.loc[df.enhancerOverlap == 1, 'alignConsensNumBases'], color='r',  alpha=1,  marker="+", s=50, label="Enhancer Overlap (n= "+str(df.loc[df.enhancerOverlap == 1, 'genomicNumBases'].shape[0])+")")
plt.title("QC: Aligment as Funciton of Enhancer Overlap")
plt.xlabel("Number of Bases in Genomic "+element) 
plt.ylabel("Number of Bases in Aligned Consensus "+element) 
plt.legend(loc=2)
plt.plot([0, max(max(df.genomicNumBases), max(df.alignConsensNumBases))], [0, max(max(df.genomicNumBases), max(df.alignConsensNumBases))], '--k') #diagonal line 
plt.show()