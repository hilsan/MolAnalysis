import pandas as pd
import os, sys
ROOT_DIR = os.path.dirname(os.path.abspath('/scratch/work/sandsth2/SCRIPTS/MolAnalysis'))
sys.path.append(ROOT_DIR)
from MolAnalysis.fingerprints import *
import joblib as jl

#Get paths to dumps of datasets to be loaded. 
path1 = sys.argv[1]
path2 = sys.argv[2]
#Fingerprint type, must be same as column
fptype = sys.argv[3]
# name of outfile 
outfile = sys.argv[4]

#Get the corresponding datasets
df1 = jl.load(path1)
df2 = jl.load(path2)

sim_df = pd.DataFrame()

for molfp in df1[fptype]:
	 sim_df[1] = df2[fptype].apply(lambda x: DataStructs.FingerprintSimilarity(x,molfp, metric=DataStructs.TanimotoSimilarity) if(np.all(pd.notnull(x))) else np.NaN)
	 
jl.dump(sim_df, outfile)
