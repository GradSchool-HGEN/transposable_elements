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


# -----------
# functions
# ----------- 

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


    # match = "MER20"
    # fastafile = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"
    # outputfilePATH = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/hmm_align/get-TFmotifs/MER20_hg38_dfam.fa

    recordFlag = False
    outfilehandle = open(outputfilePATH,'a')
    at_least_one_match = False
    with open(fastafile, 'r') as infile:
        for line in infile:
            spline = line.split(":")[0]

            if spline[0] == '>': 
                recordFlag = False
            if spline[1:] == match:
                recordFlag = True
                at_least_one_match = True 
            if recordFlag:
                # print(line)
                outfilehandle.write(line)

    outfilehandle.close()
    print("Finished. See file: "+outputfilePATH)

#-------
# Main function run 
#-------


if __name__ == "__main__":
    main(sys.argv[1:])
    #create output folder and put contents there 

