#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: ...⇧+⌘+I

import os 
import subprocess


HG38_FASTA_FILE = "/dors/capra_lab/data/dna/hg38.fasta"
REPEATMASKER_FILE = "/dors/capra_lab/data/transposable_elements/repeatmasker/sorted_filtered_hg38-TE-repeatMasker.tsv"
FIMO_MOTIFDATABASE_FILE = "/dors/capra_lab/data/TF_motifs/meme_motif_databases.12.15/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"
FIMO_OUTPUT_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_output" 
FIMO_BACKGROUND_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_background/" 
FIMO_INPUT_FASTA_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker/fimo_input"
TEMP_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/data/fimo_repeatmasker"

element = "L1MC5a"
# initialize file names for fimo
fimo_inputfasta_file = os.path.join(FIMO_INPUT_FASTA_DIR,"{}_hg38_repeatmasker.fa".format(element))
fimo_markovbackground_file = os.path.join(FIMO_BACKGROUND_DIR,"fimo_markovBackground_{}".format(element)) 
fimo_output_dir = os.path.join(FIMO_OUTPUT_DIR, "fimo_output_{}".format(element))
TE_coordinates_bedfile = temp_TE_coords.name

subprocess.run( ("fimo --verbosity 2 --oc {} --bgfile {} {} {}".format(fimo_output_dir, fimo_markovbackground_file, FIMO_MOTIFDATABASE_FILE, fimo_inputfasta_file)).split()) 
print("finished fimo on {}. it took {} seconds\n".format(round(element, time.time() - starttime, 2)))