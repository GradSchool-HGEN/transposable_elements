#!/bin/bash
## extracts coordiantes based on BED file from .maf file then filter on species 


## Abin Abraham
## 2018-03-03 23:52:10

export PATH="$PATH:/dors/capra_lab/opt/kent_tools" 

### FILE PATHS 
TE_BED_FILE="/dors/capra_lab/abraha1/projects/transposable_elements/data/dfam/MER20coords.tsv" 
OUPUT_DIR="/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative"
INPUT_MAF_DIR="/dors/capra_lab/data/alignments/hg19/multiz46way/"
CHAIN_FILE="/dors/capra_lab/abraha1/bin/liftOver_files/hg38ToHg19.over.chain"

species_to_keep="hg19\nrheMac2\nmm9"

# CHR=${INPUT_MAF_DIR##*/}
# CHR=${CHR%.maf}


### MAIN
 
printf $species_to_keep > $OUPUT_DIR/species.lst

## convert hg38 to hg19
liftOver $TE_BED_FILE $CHAIN_FILE $OUPUT_DIR/hg19_MER20_dfam_coord.tsv $OUPUT_DIR/hg38_to_hg19_unmappedLiftOver_MER20_dfam_coord.tsv

for i in {1..22}; do 
    MAF_FILE=${INPUT_MAF_DIR}chr$i.maf

    ### extract regions in BED file 
    mafsInRegion $OUPUT_DIR/hg19_MER20_dfam_coord.tsv $OUPUT_DIR/MER20_${MAF_FILE##*/} $MAF_FILE

    ### subselect species 
    mafSpeciesSubset $OUPUT_DIR/MER20_${MAF_FILE##*/} $OUPUT_DIR/species.lst $OUPUT_DIR/MER20_${MAF_FILE##*/}_speciesFiltered.maf
    rm $OUPUT_DIR/MER20_${MAF_FILE##*/}
done 

