#!/bin/python
# This script will calculate the regulatory potential for each TE occurance in the genome.
#   Regulatory Potential = (num bases with TF motif hit) / (num bases without TF motif hit)
# 
#   Input: 
#       1) -f will take file name of modified fimo output file 
#               - see for eg: /dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output/individualTE_fimo_output/ 
#       2) -o output file absolute path 
#   Output: tsv file with TE, TEsequenceID, regulatory score, number of TF hits for given sequence
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

## LOAD DATA 
df = pd.read_csv(TE_FIMO_FILE, sep='\t', delimiter=None, header=None)
df.columns = ["TF", "TE", "TEsequenceID", "TF_start", "TF_end", "strand", "motif_score",
              "p_value", "motif_seq", "q_value", "num_bases_in_TE", "num_bases_in_TFmotif"]

uniq_TE = df.TE.unique()[0]
this_TE_num_bases = int(df.loc[:, ["num_bases_in_TE"]].values[0])
uniq_seq = df.loc[:, 'TEsequenceID'].unique()

## CALC REGULATORY POTENTIAL 
output_file = open(OUTPUT_FILE, "w")
for this_seq in uniq_seq:
        motif_hit_intervals = df.loc[df['TEsequenceID'] == this_seq, ["TF_start", "TF_end"]].values
        num_TF_hits = len(motif_hit_intervals)

        # calculate number of bases with a TF hit (account for overlap of TF hits)
        num_TF_bases = numbases_TFhits(motif_hit_intervals)

        # calc regulator potential for this TE
        regulatory_potential = num_TF_bases/this_TE_num_bases
        output_file.write("{}\t{}\t{}\t{}\n".format(uniq_TE, this_seq, regulatory_potential,num_TF_hits))

output_file.close()


print("\nFinished calculating regulatory potential. Results save in:\n {}".format(OUTPUT_FILE))
print("\nTook {} seconds to finish\n".format(round(time.time()-start)))
