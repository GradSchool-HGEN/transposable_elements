#!/bin/python
# Supporting script for TE_annotate_hmm.py          
# This script will bedtools intersect dfam list of TEs with an annotation file. Output in format for TE_annotate_hmm.py
#
#
# Abin Abraham
# created on: 2018-01-03 15:05:26

import subprocess
import os
import sys 
import argparse
import datetime


TE_BED_FILE = "/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy"

#-------
# functions
#-------
def checkInputs(anno_FILE):
    '''if annotation file has > 4 column, extra column are removed, and new annotation file is returned'''
    numOfColumns = subprocess.check_output(['awk', '{print NF; exit}', anno_FILE], universal_newlines=True).splitlines()

    if int(numOfColumns[0]) > 4:
        print("\nAnnotation file has >4 columns. Extra columns will not be included for bedtools intersect.\n")
        cmd1 = 'cut -f1,2,3,4 {}'.format(anno_FILE)
        
        tidy_anno_file = "temp_{}.tsv".format(os.path.split(anno_FILE)[1])
        
        fh = open(tidy_anno_file, "w")
        subprocess.call(cmd1, stdout=fh,shell=True)
        fh.close()

    else: 
        tidy_anno_file = anno_FILE

    return tidy_anno_file

def create_annoInter(element, anno_FILE, TEbed_FILE, annoInterOut_FILE):
    anno_FILE = checkInputs(anno_FILE)
    
    filterElement = subprocess.Popen(['grep', '-w', element, TEbed_FILE], stdout=subprocess.PIPE, universal_newlines=True)
    intersectOutput = subprocess.check_output(['bedtools', 'intersect', '-wb', '-a', anno_FILE, '-b', 'stdin', '-sorted'], stdin = filterElement.stdout, stderr=subprocess.STDOUT, universal_newlines=True)  

    #write output to file
    with open(annoInterOut_FILE,'w') as f:
        f.write(intersectOutput)
    print("File Saved in: "+annoInterOut_FILE+"\n")

    return annoInterOut_FILE

#-------
# main
#-------

if __name__ == "__main__":
    # print header

    annoInterOut_FILE = os.getcwd() +"/inter-TE-anno.tsv"


    parser = argparse.ArgumentParser(description="This script runs bedtools intersect using dfamTE and annotation file")
    parser.add_argument("element", default='None',
                    action='store', type=str,
                    help="TE of interest") 
    parser.add_argument("annoFile", default='None',
                    action='store', type=str,
                    help="annoFile with Path") 
    parser.add_argument("-o", type=str, action='store',
                    default=annoInterOut_FILE, dest='outFile', 
                    help='output file path w/filename, defaults to current directory')  

    argsIN = parser.parse_args()
    elem = argsIN.element
    anno_FILE= argsIN.annoFile
    annoInterOut_FILE = argsIN.outFile

    

    # output
    print('Running {:s} {:s}'.format(' '.join(sys.argv), str(datetime.datetime.now())[:20]))
    print()
    create_annoInter(elem, anno_FILE, TE_BED_FILE, annoInterOut_FILE)





