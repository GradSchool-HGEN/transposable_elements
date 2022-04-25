#!/bin/python
# This script will retrieve fasta sequences for specific TE from dfam
#
#
#
# Abin Abraham
# created on: 2017-12-27 13:04:57

import argparse
import datetime
import sys 
import os


# -----------
# functions
# ----------- 
def extract_sequences(element_to_match, fasta_file_path, output_file_path):
    recordFlag = False

    # outputfile to write to 
    outfilehandle = open(output_file_path,'w')

    with open(fasta_file_path, 'r') as infile:
        for line in infile:
            spline = line.split(":")[0]

            if spline[0] == '>': 
                recordFlag = False

            if spline[1:] == element_to_match:
                recordFlag = True
            
            if recordFlag:
                # print(line)
                outfilehandle.write(line)

    outfilehandle.close()


    if os.path.getsize(output_file_path) == 0: 
        os.remove(output_file_path)
        print("\nAttempted to create fasta sequence for {}; but no matches were found, fasta file was not created...".format(element_to_match))
    else: 
        print("Finished creating {}".format(output_file_path))


def main(argv):

    # -----------
    # Arguments 
    # ----------- 
    parser = argparse.ArgumentParser(description="extract fasta sequence from dfam fasta file for a given element")
    parser.add_argument("TE_ofInterest", default='None',
                    action='store', type=str,
                    help="Select a transposable element") 
    parser.add_argument("inputFastaFile", 
                    action='store', type=str, 
                    help="Path for input fasta file")
    parser.add_argument("outputFilePath", 
                    action='store', type=str, 
                    help="Path for outputfile")
    
    
    argsIN = parser.parse_args()
    match = argsIN.TE_ofInterest
    outputfilePATH = argsIN.outputFilePath
    fastafile = argsIN.inputFastaFile

    extract_sequences(match, fastafile, outputfilePATH)
    # match = "MER20"
    # fastafile = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"
    # outputfilePATH = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align/get-TFmotifs/MER20_hg38_dfam.fa


#-------
# Main function run 
#-------


if __name__ == "__main__":
    main(sys.argv[1:])
    #create output folder and put contents there 

