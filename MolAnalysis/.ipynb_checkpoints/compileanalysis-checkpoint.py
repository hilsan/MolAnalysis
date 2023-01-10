import seaborn as sns
import matplotlib.pyplot as plt
from  rdkit import Chem
from rdkit.Chem import (Descriptors, AllChem, rdchem, rdmolops)
import pandas as pd
import os
from collections import defaultdict
import numpy as np

def CompileAnalysis(dataframe, foldername):
        
        if os.path.isdir(foldername):
            pass
        else:
            os.mkdir(foldername)
        stats=pd.concat([dataframe.mean(numeric_only=True),dataframe.std(numeric_only=True)], axis=1)
        stats.rename(columns={0:'mean',1:'std'},inplace=True)
        stats.to_csv(foldername+'/stats.csv',index_label="Attribute")
        plotMw(dataframe,foldername)
        plotNumberOfHeavyAtoms(dataframe,foldername)
        plotAtomicRatios(dataframe, foldername)
        print('Finished analysis in '+ foldername)

def annotateMeanStd(data,stats,string, **kws):
    mean = stats['mean'].loc[string]
    std = stats['std'].loc[string]
    ax = plt.gca()
    ax.text(.1, .6, "Mean: %.1f \n std: %.1f" % (mean,std), transform=ax.transAxes)

def plotMw(dataframe,foldername):
    sns.set_palette(sns.color_palette("bright"))
    fig, ax = plt.subplots(figsize=(7.0,4.0));
    tmp = pd.DataFrame()
    tmp['Mw'] = dataframe[dataframe.columns[dataframe.columns.get_level_values(1) == 'Mw']]
    mwplot = sns.histplot(tmp, x='Mw', bins=50, stat='percent')
    ax.set_ylabel('Percent of dataset')
    ax.set_xlabel('Molar weight [Da]')
    plt.tight_layout()
    plt.savefig(foldername+'/'+'Mw'+'_dist.png',dpi=320)
    plt.close()
        
def plotNumberOfHeavyAtoms(dataframe,foldername):
    sns.color_palette("bright")
    fig, ax = plt.subplots(figsize=(7.0,4.5));
    nAtoms = sns.histplot(dataframe, x='# heavy atoms', bins=25, discrete=True)
    ax.set_ylabel('Number of molecules')
    ax.set_xlabel('Number of heavy atoms')
    plt.tight_layout()
    plt.savefig(foldername+'/'+'heavyAtoms'+'_dist.png',dpi=320)
    plt.close()
    
    tmp = dataframe[dataframe.columns[dataframe.columns.get_level_values(0) == 'Atoms'].intersection(dataframe.columns[dataframe.columns.get_level_values(1) != 'H'])].stack()
    tmp['type'] = tmp.index.get_level_values(1)
    g = sns.catplot(tmp, x='type',y="Atoms", kind='bar')
    g.set(xlabel="Element", ylabel="Count")
    plt.tight_layout()
    plt.savefig(foldername+'/'+'heavyAtoms'+'val_dist.png',dpi=320)
    plt.close()
    
def plotAtomicRatios(dataframe, foldername):
    tmp = dataframe[dataframe.columns[dataframe.columns.get_level_values(0) == 'Atom ratios'].intersection(dataframe.columns[dataframe.columns.get_level_values(1) != 'H'])].stack()
    tmp['type'] = tmp.index.get_level_values(1)
    g = sns.catplot(tmp, x='type',y="Atom ratios", kind='box')
    g.set(xlabel="Ratio", ylabel="")
    g.set_titles(row_template="{row_name}")
    plt.tight_layout()
    plt.savefig(foldername+'/atomic_ratios.png')
    plt.close()

