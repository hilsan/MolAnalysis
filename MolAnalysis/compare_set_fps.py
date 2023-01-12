"Function to compare fingerprints from two different datasets"
import os
import sys
import joblib as jl
import numpy as np
import pandas as pd
ROOT_DIR = os.path.dirname(os.path.abspath('/scratch/work/sandsth2/SCRIPTS/MolAnalysis'))
sys.path.append(ROOT_DIR)
from MolAnalysis.fingerprints import *


#Get paths to dumps of datasets to be loaded.

PATH_1 = sys.argv[1]
PATH_2 = sys.argv[2]
#Fingerprint type, must be same as column
fptype = sys.argv[3]
# name of outfile
outfile = sys.argv[4]

#Get the corresponding datasets
df1 = jl.load(PATH_1)
df2 = jl.load(PATH_2)

sim_df = pd.DataFrame()

for molfp in df1[fptype]:
    sim_df[1] = df2[fptype].apply(lambda x:
    	DataStructs.FingerprintSimilarity(x, molfp, metric=DataStructs.TanimotoSimilarity)
    	if(np.all(pd.notnull(x))) else np.NaN)

jl.dump(sim_df, outfile)
