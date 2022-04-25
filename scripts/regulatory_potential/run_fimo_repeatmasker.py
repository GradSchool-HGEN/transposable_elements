#!/bin/python
# This script will run fimo on each TE in repeatmasker (after filtering certain elements). 
#   > inputs: 
#       - 
#       -  
#
#   > depends: 
#       - fimo utility from meme suite 
#       - 
#
#   > outputs: 
#       - fimo_<TEName>.txt: contains fimo-output with qvalues 
#               - columns are: TF, TE, TEsequenceID, TF_start(1based,closed), TF_end(1based_closed), strand, motif_score, p_value, motif_seq, q_value, num_bases_in_TE, num_bases_in_TFmotif
#
# abin-personal Abraham
# created on: 2018-01-16 13:00:08
# updated on: 2018-01-18 09:27:18

import subprocess 
import os
import tempfile
import time 
import pandas as pd 
import numpy as np 
from multiprocessing import Pool 
from functools import partial

starttime = time.time() 
# /Users/abin-personal-personal/Desktop/transfer/scripts/dev_calc_regPotential.py
#-------
# dependencies
#-------
# HG38_FASTA_FILE = "/Users/abin-personal/Desktop/transfer/data/hg38.fasta"
# REPEATMASKER_FILE = "/Users/abin-personal/Desktop/transfer/data/sorted_filtered_hg38-TE-repeatMasker.tsv"
# FIMO_MOTIFDATABASE_FILE = "/Users/abin-personal/Desktop/transfer/data/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
# FIMO_OUTPUT_DIR = "/Users/abin-personal/Desktop/transfer/output/fimo_out/" 
# FIMO_BACKGROUND_DIR = "/Users/abin-personal/Desktop/transfer/output/fimo_background/" 
# FIMO_INPUT_FASTA_DIR = "/Users/abin-personal/Desktop/transfer/inputs/element_fasta"
# TEMP_DIR = "/Users/abin-personal/Desktop/transfer/scripts"


HG38_FASTA_FILE = "/dors/capra_lab/data/dna/hg38.fasta"
REPEATMASKER_FILE = "/dors/capra_lab/data/transposable_elements/repeatmasker/sorted_filtered_hg38-TE-repeatMasker.tsv"
FIMO_MOTIFDATABASE_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/HOCOMOCOv10_HUMAN_mono_meme_format.meme_filtered"
FIMO_OUTPUT_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output" 
FIMO_BACKGROUND_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_background" 
FIMO_INPUT_FASTA_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_input"
TEMP_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/temp"

#-------
# functions
#-------
def get_fasta_from_bed(fimo_inputfasta_file, HG38_FASTA_FILE, TE_coordinates_bedfile):
    subprocess.run(("bedtools getfasta -fo {} -fi {} -bed {}".format(fimo_inputfasta_file, HG38_FASTA_FILE, TE_coordinates_bedfile)).split())

def get_markovbackground(fimo_inputfasta_file, fimo_markovbackground_file):
    subprocess.run( ("fasta-get-markov {} > {}".format(fimo_inputfasta_file, fimo_markovbackground_file)), shell=True)

def check_fimo_parameters(fimo_inputfasta_file, fimo_markovbackground_file, TE_coordinates_bedfile):
    if os.path.isfile(fimo_inputfasta_file) == False: 
         get_fasta_from_bed(fimo_inputfasta_file, HG38_FASTA_FILE, TE_coordinates_bedfile)

    if os.path.isfile(fimo_markovbackground_file) == False: 
        get_markovbackground(fimo_inputfasta_file, fimo_markovbackground_file)

    # if os.path.isdir(FIMO_OUTPUT_DIR) == False: 
    #     os.mkdir(FIMO_OUTPUT_DIR)

def run_fimo(df, element):

    ## filter, create TE coord file 
    one_TE_coords_df = df.loc[df['TE'] == element, ['chr','start','end']]
    temp_TE_coords_file = tempfile.NamedTemporaryFile(mode='w', newline='\n', dir=TEMP_DIR)
    one_TE_coords_df.to_csv(temp_TE_coords_file.name, sep='\t', header=False, index=False)

    ## create file names for fimo
    fimo_inputfasta_file = os.path.join(FIMO_INPUT_FASTA_DIR,"{}_hg38_repeatmasker.fa".format(element))
    fimo_markovbackground_file = os.path.join(FIMO_BACKGROUND_DIR,"fimo_markovBackground_{}".format(element)) 
    TE_coordinates_bedfile = temp_TE_coords_file.name
    
    check_fimo_parameters(fimo_inputfasta_file, fimo_markovbackground_file, TE_coordinates_bedfile) 

    ## run fimo
    print("*'{}': processing {} and about to run fimo".format(os.path.basename(__file__), element))
    fimo_start_time = time.time()
    try: 
        # this_fimo_output = subprocess.check_output( ("fimo --text --skip-matched-sequence --motif AHR_HUMAN.H10MO.B --verbosity 2 --bgfile {} {} {}".format(fimo_markovbackground_file, FIMO_MOTIFDATABASE_FILE, fimo_inputfasta_file)), shell = True,universal_newlines=True) 
        this_fimo_output = subprocess.check_output( ("fimo --text --verbosity 1 --bgfile {} {} {}".format(fimo_markovbackground_file, FIMO_MOTIFDATABASE_FILE, fimo_inputfasta_file)), shell = True,universal_newlines=True) 
    except subprocess.CalledProcessError as er: 
        return print(er)

    temp_file = open(os.path.join(TEMP_DIR,"temp_fimoout_{}".format(element)), mode='w')
    temp_file.writelines(this_fimo_output)
    temp_file.close()
    
    ## generate q-values
    try: 
        qval_output = subprocess.check_output("qvalue --header 1 --column 7 --append {}".format(os.path.join(TEMP_DIR,"temp_fimoout_{}".format(element))), stderr=subprocess.STDOUT, shell=True, universal_newlines=True )
    except subprocess.CalledProcessError as er2: 
        return print(er2) 

    qval_list = qval_output.splitlines() 

    fimo_complete_time = time.time()
    print("*'{}': finished fimo w/ qvalues on {}. it took {} seconds\n".format(os.path.basename(__file__), element, round(fimo_complete_time - fimo_start_time, 2)))
    
    ## write fimo results w/ qvalue to file 
    ## append number of bases in TF hit and TE sequence and append TE name 
    ## final output file has following columns: TF, TE, TEsequenceID, TF_start(1based,closed), TF_end(1based_closed), strand, motif_score, p_value, motif_seq, q_value, num_bases_in_TE, num_bases_in_TFmotif
    with open(os.path.join(FIMO_OUTPUT_DIR, "fimo_{}.txt".format(element)), 'w') as outfile:
        for i,this_row in enumerate(qval_list): 
            if i < 2: 
                continue

            genome_seq_coord = this_row.split()
            start, end = genome_seq_coord[1].split(":")[1].split('-')
            TFstart, TFend = genome_seq_coord[2:4]
            genome_seq_coord.insert(1, element)
            genome_seq_coord.append(str(int(end)-int(start)))
            genome_seq_coord.append(str(int(TFend) - int(TFstart) + 1))
            if float(genome_seq_coord[7]) < 0.10:
                outfile.write('\t'.join(genome_seq_coord)+'\n')

    
    print("*'{}': finished writing fimo-output-with-qVal to {}. Total time was {} seconds\n".format(os.path.basename(__file__), os.path.join(FIMO_OUTPUT_DIR, "fimo_{}.txt".format(element)), round(time.time() - starttime, 2)))
#-------
# main
#-------

df = pd.read_csv(REPEATMASKER_FILE, sep='\t', header=None)
df.columns = ["chr","start","end","TE","Family"]
uniq_te=df.loc[~df.duplicated('TE', keep = 'first'),'TE']


# element = "L1MC5a"
#run_fimo(df, element)

pool = Pool()
partial_run_fimo = partial(run_fimo, df)


exp_sum_list = pool.map(partial_run_fimo, uniq_te)
pool.close()
pool.join()


