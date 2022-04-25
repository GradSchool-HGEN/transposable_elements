#!/bin/python
# This script will examine the alignment of TE elements in Dfam. 
#   Plots num of bases in aligned consensus vs. aligned genomic locus for a given TE 
#   
#
# Abin Abraham
# created on: 2017-12-16 11:34:06


from Bio import SeqIO
import subprocess
import pandas as pd 
import numpy as np 
import editdistance as ed 
import matplotlib.pyplot as plt 
import matplotlib.path as path
import matplotlib.patches as patches
# from multiprocessing import Pool 

# -----------
# Functions
# ----------- 


def TE_rename(id):
    newName = id.split(":")[2:4]
    newName = '-'.join(newName)[:-3]
    return newName

def consensus_renameKey(id):
    newName = id.split("-")[0:-1]
    newName = '-'.join(newName)
    return newName

def getConsSeq(element):
    consensus_FILE = "/dors/capra_lab/data/transposable_elements/dfam/Dfam.cons.fa"
    consSeq = SeqIO.index(consensus_FILE, "fasta", key_function = consensus_renameKey)     
    consString = consSeq[element].seq
    return consString

# -----------
# main
# ----------- 

element="MER20"

# load dfam TE coordinates & filter for element & load in to pandas
TEbedfile_PATH = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"
filterForElement = subprocess.check_output(['grep', '-w', element, TEbedfile_PATH], universal_newlines=True)
out_data =[x.split('\t') for x in filterForElement.split('\n')]
col_name = ["chr", "start", "end", "TE", "model_start","model_end","model_length"]
df = pd.DataFrame(data=out_data[:-1], columns=col_name)
df = df.apply(pd.to_numeric, errors='ignore')

#calcualte number of bases in aligned-consensus & aligned genomic hit 
genomicNumBases = df.end - df.start 
alignConsensNumBases = df.model_end - df.model_start
deviationFromConsensus = abs(genomicNumBases-alignConsensNumBases)
per_deviationFromConsensus = abs(genomicNumBases-alignConsensNumBases)/df.model_length

#plot 
plt.scatter(genomicNumBases, alignConsensNumBases,alpha=0.6,  marker=".", s =2, label="one "+element+" locus")
plt.title("Comparing Dfam Alignment of "+element) 
plt.xlabel("Number of Bases in Genomic "+element) 
plt.ylabel("Number of Bases in Aligned "+element) 
plt.legend(loc=2)
plt.plot([0, max(max(genomicNumBases), max(alignConsensNumBases))], [0, max(max(genomicNumBases), max(alignConsensNumBases))], '--r') #diagonal line 

plt.show()

#plot histogram of deviation from consensus 
fig, ax = plt.subplots()
n, bins = np.histogram(per_deviationFromConsensus, 100)
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n
XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T
barpath = path.Path.make_compound_path_from_polys(XY)
patch = patches.PathPatch(barpath)
ax.add_patch(patch)
ax.set_xlim(left[0], right[-1])
ax.set_ylim(bottom.min(), top.max())


plt.title("Deviation of alignment of " + element+ "to Consensus") 
plt.xlabel("Base Difference between Consensus and Genomic Hit Divided by " +element+ " length")
plt.ylabel("Count") 
plt.show()








# code to calculate hamming/Levenstein Distance below

# TESeq_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"
# TESeq = SeqIO.index(TESeq_FILE, "fasta", key_function=TE_rename)
# EDVal = np.empty(len(out_data))


# for i,row in df.iterrows(): 
#     consStart = int(row.model_start)
#     consEnd =  int(row.model_end)
#     elemString = '-'.join(df.iloc[i, [0,1,2]].values)
#     consString = str(getConsSeq(element)[consStart:consEnd])
#     EDVal[i] = ed.eval(elemString, consString)
#     print(i)


# EDVal = list()
# for i,x in enumerate(out_data):
#     print(i)
#     consStart = int(out_data[i][4])
#     consEnd =  int(out_data[i][5])
#     elemString = str(TESeq['-'.join(out_data[i][0:3])].seq)
#     consString = str(getConsSeq(element)[consStart:consEnd])

#     EDVal.append(ed.eval(elemString, consString))

# if False: 

#     import matplotlib.mlab as mlab
#     import matplotlib.pyplot as plt

#     mu, sigma = 100, 15
#     x = EDVal

#     # the histogram of the data
#     n, bins, patches = plt.hist(x, normed=1, facecolor='green', alpha=0.75)

#     # add a 'best fit' line
#     y = mlab.normpdf( bins, mu, sigma)
#     l = plt.plot(bins, y, 'r--', linewidth=1)

#     plt.xlabel('Smarts')
#     plt.ylabel('Probability')
#     plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#     plt.axis([40, 160, 0, 0.03])
#     plt.grid(True)

#     plt.show()