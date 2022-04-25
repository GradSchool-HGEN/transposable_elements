#!/bin/python
# This script will calculate the regulatory potential for each TE for each TF on each occurance in the genome.
#   Regulatory Potential = (num bases with TF motif hit) / (num bases without TF motif hit)
#
#   Output: tsv file with TE, TEsequenceID, TFmotif, regulatory score, number of TF hits
#
#
# Abin Abraham
# created on: 2018-01-27 10:05:07

import os
import time
import argparse
import pandas as pd
import numpy as np

start = time.time()

# note: TF motif coordinates are 1based closed
ROOTPATH = "/Users/abin-personal/Desktop/transfer/data/"
TE_FIMO_FILE = "fimo_Alu.txt"
OUTPUT_FILE = "output_regPotential_byTF.tsv"


# -----------
# functions
# -----------
def numbases_TFhits(hits_interval):
    expand_indicies = [np.arange(thisTF_interval[0], thisTF_interval[1]+1)
                       for thisTF_interval in hits_interval]
    unwrap_indicies = [
        item for thisInterval in expand_indicies for item in thisInterval]
    TF_base_hit_indicies = set(unwrap_indicies)
    num_TF_bases = len(TF_base_hit_indicies)
    return num_TF_bases


# -----------
# main
# -----------

## PARSE ARGS, obtain input file path
parser = argparse.ArgumentParser(description="calculate ratio of # of bases with and without TF motif hit")
parser.add_argument("-f",
                    dest = 'file_path', action='store', type=str,
                    help="absolute path to input file") 
parser.add_argument("-o",
                    dest = 'output_file_path', action='store', type=str,
                    help="absolute path to output file") 

argsIN = parser.parse_args()
TE_FIMO_FILE = argsIN.file_path
OUTPUT_FILE = argsIN.output_file_path

print("\nRunning calc_reg_potential.py on input file: {}".format(TE_FIMO_FILE))

# LOAD DATA
df = pd.read_csv(os.path.join(ROOTPATH, TE_FIMO_FILE),
                 sep='\t', delimiter=None, header=None)
df.columns = ["TF", "TE", "TEsequenceID", "TF_start", "TF_end", "strand", "motif_score",
              "p_value", "motif_seq", "q_value", "num_bases_in_TE", "num_bases_in_TFmotif"]

uniq_TE = df.TE.unique()[0]
this_TE_num_bases = int(df.loc[:, ["num_bases_in_TE"]].values[0])
uniq_seq = df.loc[:, 'TEsequenceID'].unique()


output_file = open(os.path.join(ROOTPATH, OUTPUT_FILE), "w")

for this_seq in uniq_seq:
    uniq_TFs = df.loc[df['TEsequenceID'] == this_seq, 'TF'].unique()

    for this_TF in uniq_TFs:
        this_TE_num_bases = int(df.loc[df['TEsequenceID'] == this_seq, ["num_bases_in_TE"]].values[0])
        motif_hit_intervals = df.loc[(df['TEsequenceID'] == this_seq) & (df['TF'] == this_TF), ["TF_start", "TF_end"]].values
        num_TF = len(motif_hit_intervals)

        # calculate number of bases with a TF hit (account for overlap of TF hits)
        num_TF_bases = numbases_TFhits(motif_hit_intervals)
        # calc regulator potential for this TE
        regulatory_potential = num_TF_bases/this_TE_num_bases
        output_file.write("{}\t{}\t{}\t{}\t{}\n".format(uniq_TE, this_seq, this_TF, regulatory_potential, num_TF))

output_file.close()

print("\nFinished calculating regulatory potential **by TF**. Results save in: {}".format(OUTPUT_FILE))
print("Took {} seconds to finish\n".format(round(time.time()-start)))
