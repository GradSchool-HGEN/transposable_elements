#!/bin/bash
# this script will use the fasta file for a TE of interest and then run fimo on it 
# required inputfiles 
#   > dfam fasta file IF fasta file for given TE (element) has to be created
#   > element fasta file: fasta sequences of all sequences to be searched by fimo 

#***CHANGE THIS AS FIT ****
element="MER20" 

#dependencies
dfamFasta_file="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hg38_dfam.nrph.hits_tidy.fasta"
motifFile="/dors/capra_lab/data/TF_motifs/meme_motif_databases.12.15/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"

#inputfiles required for running fimo
elemFastaFile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/element_fasta/${element}_hg38_dfam.fa"
markovOutput="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/fimo-dfam/markov_${element}"
consensusFastaFile="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/hmmAlign/align-input_${element}/consensus_${element}.fa"
#outputfiles
fimoOutput_dir="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/fimo-dfam/fimo-out-${element}"
fimoOutputCons_dir="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/fimo-dfam/fimo-out-consensusSeq_${element}"


##
#check if fasta sequence already exists for given element 
if [ -f $elemFastaFile ]; then
   echo "Fasta Sequence File already exists...\n"
else
   echo "File $elemFastaFile does not exist, creating it..."
   python getFasta_dfamTE.py $element $dfamFasta_file $elemFastaFile
fi

#check if fasta sequence already exists for given element 
if [ -f $consensusFastaFile ]; then
   echo "Consensus Fasta Found...\n"
else
   echo "File $elemFastaFile does not exist, creating it..."
   exit -1
fi

#check if markov background file already exists for given element 
if [ -f $markovOutput ]; then
   echo "Markov Background File already exists..."
   echo ""
else
   echo "File $markovOutput does not exist, creating it..."
   echo ""
   fasta-get-markov $elemFastaFile > $markovOutput
fi


# fimo --o $fimoOutput_dir -bgfile $markovOutput $motifFile $elemFastaFile
fimo --o $fimoOutputCons_dir -bgfile $markovOutput $motifFile $consensusFastaFile
# fimo --text -bgfile $markovOutput $motifFile $elemFastaFile > $fimoOutput_dir.txt

## 
# format fimo output

cd $fimoOutput_dir

# awk -v OFS='\t' 'NR>1{split($1,TFf,"_"); split($1,Tfq,"."); split($2,seqName,"::"); split(seqName[2],chr,":"); split(chr[2],rcoord, "\(" ); split(rcoord[1],startEnd, "-"); if((Tfq[3] =="A" || Tfq[3] =="B" || Tfq[3]=="C") && ($8<0.10))print chr[1],(startEnd[1]+$3-1),(startEnd[2]+$4),$6,  TFf[1],$2,$3,$4,$5,$7,$8,$9,Tfq[3] }' fimo.txt > formatted_fimo_out_${element}.tsv
# sort -k1,1 -k2,2n formatted_fimo_out_${element}.tsv -o formatted_fimo_out_${element}.tsv


awk -v OFS='\t' 'NR>1{split($1,TFf,"_"); split($1,Tfq,"."); split($2,seqName,"::"); split(seqName[2],chr,":"); split(chr[2],rcoord, "\(" ); split(rcoord[1],startEnd, "-"); if((Tfq[3] =="A" || Tfq[3] =="B" || Tfq[3]=="C") && ($8<0.10))print chr[1],(startEnd[1]+$3-1),(startEnd[2]+$4),$6,  TFf[1],$2,$3,$4,$5,$7,$8,$9,Tfq[3] }' fimo.txt > formatted_fimo_out_${element}.tsv
sort -k1,1 -k2,2n formatted_fimo_out_${element}.tsv -o formatted_fimo_out_${element}.tsv

#formatted file columns 
# chr, start_relative, end_relative, TF, score SeqName, start, stop, strand, pvalue, qvalue,matched, TF quality