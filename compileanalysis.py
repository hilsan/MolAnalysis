"Function to gather descriptor statistics and plot them"
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def compile_analysis(dataframe, foldername):
    "Function to gather descriptor statistics and plot them and put them in foldername"
    if os.path.isdir(foldername):
        pass
    else:
        os.mkdir(foldername)
    stats=pd.concat([dataframe.droplevel(0,axis=1).mean(numeric_only=True),dataframe.droplevel(0,axis=1).std(numeric_only=True)], axis=1)
    stats.rename(columns={0:'mean',1:'std'},inplace=True)
    stats.to_csv(foldername+'/stats.csv',index_label="Attribute")

def plot_mw(dataframe, label, foldername):
    "Plot the distribution of the molar weight"
    sns.set_palette(sns.color_palette("bright"))
    cm=2.54
    fig, ax1 = plt.subplots(figsize=(6, 4), dpi=300, constrained_layout=True)
    dataframe['Mw'] = dataframe['Mw'].astype('float32')
    sns.histplot(dataframe, x='Mw', bins=40, ax=ax1, stat='percent')
    plt.ylabel('Percent')
    plt.xlabel('Molar weight [Da]')
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    plt.tight_layout()
    plt.savefig(foldername+'/'+label+'_Mw_dist.png',dpi=320)
    plt.close()

def plot_number_of_heavy_atoms(dataframe, label, foldername):
    "Plot the distribution of the heavy atoms"
    sns.color_palette("bright")
    dataframe['# heavy atoms'] = dataframe['# heavy atoms'].astype('float32')
    fig, ax1 = plt.subplots(figsize=(6, 4), dpi=300, constrained_layout=True)
    sns.histplot(dataframe, x='# heavy atoms', ax=ax1, bins=25, discrete=True)
    plt.ylabel('Number of molecules')
    plt.xlabel('Number of heavy atoms')
    ax1.set_xticks(np.arange(0, np.max(dataframe['# heavy atoms'])+1, 10))
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')

    plt.tight_layout()
    plt.savefig(foldername+'/'+label+'_heavyAtoms_dist.png',dpi=320)
    plt.close()
    
def plot_atom_distribution(dataframe, label, foldername):
    fig, ax1 = plt.subplots(figsize=(6, 4), dpi=300, constrained_layout=True)
    element_plot = sns.catplot(dataframe, x='type',y=0, kind='bar')
    element_plot.set(xlabel="Element", ylabel="Count")
    plt.tight_layout()
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    plt.savefig(foldername+'/'+label+'_heavyAtoms_val_dist.png',dpi=320)
    plt.close()

def arange_atomic_distribution(dataframe):
    df2 = dataframe.filter(regex='Atoms').copy(deep=True)
    df2.columns = df2.iloc[0]
    df2.drop(0, inplace=True)
    tmp = pd.concat([df2[df2.columns[df2.columns.isin(['C','N',
                                   'O'])]].astype('float32'), pd.DataFrame(df2[df2.columns[~df2.columns.isin(['C','N',
                                   'O','H'])]].astype('float32').sum(axis=1)).rename({0:'other'},axis='columns')],axis=0).astype('float32').stack()
    tmp.index = tmp.index.droplevel()
    tmp.index.rename('type', inplace=True)
    return tmp
    
def arange_atomic_ratios(dataframe, label):
    dr = dataframe[['O:C','H:C', 'N:C']].astype('float32').stack()
    dr.index = dr.index.droplevel(None)
    dr.index.rename('Type', inplace=True)
    a_r = pd.concat([pd.DataFrame(dr.reset_index()).rename({0:'Ratio'}, axis='columns'), 
               pd.DataFrame(len(dr)*[label]).rename({0:'label'}, axis='columns')], axis=1)
    return a_r
    
    
def plot_atomic_ratios(dataframe, foldername):
    fig, ax1 = plt.subplots(figsize=(6, 4), dpi=300, constrained_layout=True)
    plt.rc('font', size=12)
    sns.set_palette(sns.color_palette("bright"))
    sns.pointplot(data=dataframe,  ax=ax1, errwidth=1,
                  x="label", y="Ratio",
                  hue="Type",  errorbar="sd",
                  capsize=0.05, join=False)
    ax1.set_ylim((-0.2,3))
    ax1.legend(loc='center left', 
               bbox_to_anchor=(1, 0.5),
               prop={'size': 12},
               fancybox=True,
               framealpha=0)
    plt.setp(ax1.get_xticklabels(), rotation=30, horizontalalignment='right')
    plt.savefig(foldername+'/atom_ratios',transparent=True)
