#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2018-02-28 09:50:09

import os 
import pickle
import pandas as pd
import numpy as np 

WORKING_DIR = "/dors/capra_lab/abraha1/projects/transposable_elements/scripts/tf_dynamics/MER20/mapped_TF"
this_mapped_data = "MER20_CTCF_mappedData_mult_2018-02-26.npy"

mapped_dict = np.load(os.path.join(WORKING_DIR,this_mapped_data)).item()

store_hits = np.zeros(len(list(mapped_dict.values())[0]))
for this_seq, this_array in mapped_dict.items(): 
    # print(this_seq)
    print(this_array)
    np.isnan()





# CTCF 