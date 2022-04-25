#!/bin/bash
## extracts coordiantes from multiz64way based on input BED file from; then filters on species 


## Abin Abraham
# 2018-03-23 16:58:00

export PATH="$PATH:/dors/capra_lab/opt/kent_tools" 


######
### TO DO 
### modify TE_BED_FILE to include all fiel and then continue modifying scripts 

### FILE PATHS 
OUPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative"
INPUT_MAF_DIR="/dors/capra_lab/data/alignments/hg19/multiz46way/"
CHAIN_FILE="/dors/capra_lab/abraha1/bin/liftOver_files/hg38ToHg19.over.chain"

TE_BED_FILE=(AluJo_hg19.out  AluSg_hg19.out  AluSq_hg19.out  AluYa5_hg19.out  AluY_hg19.out AluJb_hg19.out  AluSc_hg19.out  AluSp_hg19.out  AluSx_hg19.out  AluYb8_hg19.out)
TE_BED_PATH="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/alu_comparative/repeatmasker_Alu/"

species_to_keep="hg19\nrheMac2\nmm9"

# CHR=${INPUT_MAF_DIR##*/}
# CHR=${CHR%.maf}


### MAIN
 
printf $species_to_keep > $OUPUT_DIR/species.lst

# ## convert hg38 to hg19
# liftOver $TE_BED_FILE $CHAIN_FILE $OUPUT_DIR/hg19_MER20_dfam_coord.tsv $OUPUT_DIR/hg38_to_hg19_unmappedLiftOver_MER20_dfam_coord.tsv

THIS_TE_FILE=${TE_BED_FILE[1]}
THIS_TE_FILE_PATH=${TE_BED_PATH}${THIS_TE_FILE}


for THIS_TE_FILE in ${TE_BED_FILE[@]}
do
THIS_TE_FILE_PATH=${TE_BED_PATH}${THIS_TE_FILE}

echo TE_FILE: $THIS_TE_FILE
    # for i in {1..22}; do 
    for i in {1..22}; do 
        MAF_FILE=${INPUT_MAF_DIR}chr$i.maf

        echo CHR: $MAF_FILE
        ### extract regions in BED file 
            # mafsInRegion regions.bed out.maf|outDir in.maf(s)
        mafsInRegion $THIS_TE_FILE_PATH $OUPUT_DIR/${THIS_TE_FILE}_${MAF_FILE##*/} $MAF_FILE

        ### subselect species 
        mafSpeciesSubset $OUPUT_DIR/${THIS_TE_FILE}_${MAF_FILE##*/} $OUPUT_DIR/species.lst $OUPUT_DIR/${THIS_TE_FILE}_${MAF_FILE##*/}_speciesFiltered.maf
        
        rm $OUPUT_DIR/${THIS_TE_FILE}_${MAF_FILE##*/}
        echo "check" $OUPUT_DIR/${THIS_TE_FILE}_${MAF_FILE##*/}_speciesFiltered.maf
    done 


done

