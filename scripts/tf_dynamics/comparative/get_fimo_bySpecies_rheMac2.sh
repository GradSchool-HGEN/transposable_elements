#!/bin/bash
#  This script runs fimo on MER20s across human, and homologs sequences of Rhesus and Mouse. 
# 
# 
# 
# required inputfiles 
#   > dfam fasta file IF fasta file for given TE (element) has to be created
#   > element fasta file: fasta sequences of all sequences to be searched by fimo 


##
# FILE PATHS -- USER MUST MODIFY 
##

hg19_FASTA_FILE="/dors/capra_lab/abraha1/data/dna/hg19.fa"
mm9_FASTA_FILE="/dors/capra_lab/data/dna/mm9/mm9.fa"
rheMac2_FASTA_FILE="/dors/capra_lab/data/dna/rheMac2/rheMac2.fa"

hg19_MER20_FASTA_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/fasta_files/hg19_MER20.fa"
mm9_MER20_FASTA_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/fasta_files/mm9_MER20.fa"
rheMac2_MER20_FASTA_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/fasta_files/rheMac2_MER20.fa"

hg19_MER20_COORD_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/hg19_MER20_coords.tsv"
mm9_MER20_COORD_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/mm9_MER20_coords.tsv"
rheMac2_MER20_COORD_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/rheMac2_MER20_coords.tsv"

MEME_MOTIF_DATABASE="/dors/capra_lab/data/TF_motifs/meme_motif_databases.12.15/JASPAR/JASPAR_CORE_2016_vertebrates.meme"

hg19_MARKOV_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/markov/hg19_markov"
mm9_MARKOV_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/markov/mm9_markov"
rheMac2_MARKOV_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/markov/rheMac2_markov"

hg19_FIMO_OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/hg19_fimo_out"
mm9_FIMO_OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/mm9_fimo_out"
rheMac2_FIMO_OUTPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/intersect_fimo/rheMac2_fimo_out"


##
# GET FASTA 
##

# if [ -f $hg19_MER20_FASTA_FILE ]; then
#    echo "Fasta Sequence File already exists...\n"
# else
#    echo "File $hg19_MER20_FASTA_FILE does not exist, creating it..."
#    bedtools getfasta -fi $hg19_FASTA_FILE -bed $hg19_MER20_COORD_FILE -fo $hg19_MER20_FASTA_FILE
# fi

# if [ -f $mm9_MER20_FASTA_FILE ]; then
#    echo "Fasta Sequence File already exists...\n"
# else
#    echo "File $mm9_MER20_FASTA_FILE does not exist, creating it..."
#    bedtools getfasta -fi $mm9_FASTA_FILE -bed $mm9_MER20_COORD_FILE -fo $mm9_MER20_FASTA_FILE
# fi

if [ -f $rheMac2_MER20_FASTA_FILE ]; then
   echo "Fasta Sequence File already exists...\n"
else
   echo "File $rheMac2_MER20_FASTA_FILE does not exist, creating it..."
   bedtools getfasta -fi $rheMac2_FASTA_FILE -bed $rheMac2_MER20_COORD_FILE -fo $rheMac2_MER20_FASTA_FILE
fi




##
# GET MARKOV BACKGROUND 
##

# if [ -f $hg19_MARKOV_FILE ]; then
#    echo "Markov Background File already exists...\n"
# else
#    echo "File $hg19_MARKOV_FILE does not exist, creating it...\n"
#    fasta-get-markov $hg19_FASTA_FILE > $hg19_MARKOV_FILE
# fi

# if [ -f $mm9_MARKOV_FILE ]; then
#    echo "Markov Background File already exists...\n"
# else
#    echo "File $mm9_MARKOV_FILE does not exist, creating it...\n"
#    fasta-get-markov $mm9_FASTA_FILE > $mm9_MARKOV_FILE
# fi

if [ -f $rheMac2_MARKOV_FILE ]; then
   echo "Markov Background File already exists...\n"
else
   echo "File $rheMac2_MARKOV_FILE does not exist, creating it...\n"
   fasta-get-markov $rheMac2_FASTA_FILE > $rheMac2_MARKOV_FILE
fi


##
# RUN FIMO 
##

# fimo --o $hg19_FIMO_OUTPUT_DIR -bgfile $hg19_MARKOV_FILE $MEME_MOTIF_DATABASE $hg19_MER20_FASTA_FILE
# fimo --o $mm9_FIMO_OUTPUT_DIR -bgfile $mm9_MARKOV_FILE $MEME_MOTIF_DATABASE $mm9_MER20_FASTA_FILE
fimo --o $rheMac2_FIMO_OUTPUT_DIR -bgfile $rheMac2_MARKOV_FILE $MEME_MOTIF_DATABASE $rheMac2_MER20_FASTA_FILE






##
#check if fasta sequence already exists for given element 
# if [ -f $elemFastaFile ]; then
#    echo "Fasta Sequence File already exists...\n"
# else
#    echo "File $elemFastaFile does not exist, creating it..."
#    python getFasta_dfamTE.py $element $dfamFasta_file $elemFastaFile
# fi

# #check if fasta sequence already exists for given element 
# if [ -f $consensusFastaFile ]; then
#    echo "Consensus Fasta Found...\n"
# else
#    echo "File $elemFastaFile does not exist, creating it..."
#    exit -1
# fi