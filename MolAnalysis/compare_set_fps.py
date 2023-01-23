#!/usr/bin/env python
"Function to compare fingerprints from two different datasets"
import os
import sys
import joblib as jl
import time
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import (MACCSkeys, RDKFingerprint, DataStructs)
ROOT_DIR = os.path.dirname(os.path.abspath('/scratch/work/sandsth2/SCRIPTS/MolAnalysis'))
sys.path.append(ROOT_DIR)
from MolAnalysis.fingerprints import *

datasets={'Gecko':'/scratch/work/sandsth2/DataSets/gecko/gecko_full.dump','Wang':'/scratch/work/sandsth2/DataSets/acp-17-7529-2017-supplement/wang_data.dump','QM9':"/scratch/work/sandsth2/DataSets/QM9/qm9_data.dump",'Quinones':"/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_nosulf.dump",'ExpMoNA':"/scratch/work/sandsth2/DataSets/MassBank/expmona_cleaned.dump",'nablaDFT':"/scratch/work/sandsth2/DataSets/nablaDFT/nablaDFT.dump"}
#--------Functions------------------------------------------------
def similarity(molfp,df2): 
    return df2['MACCS'].apply(lambda x:DataStructs.FingerprintSimilarity(
            x, molfp, metric=DataStructs.TanimotoSimilarity) if(np.all(pd.notnull(x))) else np.NaN)
def get_df(key): 
    name = key+'_maccs_fp.dump'
    if os.path.isfile(name):
        df = jl.load(name)
    else:
        df = jl.load(datasets[key])
        MACCSAnalysis(df)
        jl.dump(df, name,compress=3)
    return df

def get_batch_df(key): 
    name = key+'.dump'
    try:
        df = jl.load(name)
        MACCSAnalysis(df)
    except:
        print('Failed loading batch dataframe.')
    return df
#------------Set constants and variables---------------------
n_jobs = int(sys.argv[1])
print(sys.argv)
key1 = sys.argv[2]
key2 = sys.argv[3]

df1 = get_df(key1)
df2 = get_batch_df(key2)

#------------Main-------------------------
#Fingerprint type, must be same as column
fptype='MACCS'

print('Starting analysis on '+key2)
# name of outfile
if os.path.isdir('comp'):
    pass
else:
    os.mkdir('comp')

outfile = 'comp'+"/"+key2+'_'+key1+'.dump'
#-------------------------------------
t = time.time()
sim_df = pd.DataFrame()
results = Parallel(n_jobs=n_jobs)(delayed(similarity)(molfp,df2) for molfp in df1['MACCS'].values)
sim_df = pd.concat([sim_df, pd.DataFrame(results)], axis = 1)
print(time.time() - t)

jl.dump(sim_df.reset_index(), outfile,compress=4)
#-----plot----------------------------  
#ax = sim_df.stack().plot.hist()
#fig = ax.get_figure()
#fig.savefig(outfile+'.png')
#plt.close()
