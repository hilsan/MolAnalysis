"Function to gather descriptor statistics and plot them"
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def compile_analysis(dataframe, foldername):
    "Function to gather descriptor statistics and plot them and put them in foldername"
    if os.path.isdir(foldername):
        pass
    else:
        os.mkdir(foldername)
    stats=pd.concat([dataframe.mean(numeric_only=True),dataframe.std(numeric_only=True)], axis=1)
    stats.rename(columns={0:'mean',1:'std'},inplace=True)
    stats.to_csv(foldername+'/stats.csv',index_label="Attribute")
    plot_mw(dataframe,foldername)
    plot_number_of_heavy_atoms(dataframe,foldername)
    plot_atomic_ratios(dataframe, foldername)
    print('Finished analysis in '+ foldername)

def annotate_mean_std(stats,string):
    "Write mean, std to figure"
    mean = stats['mean'].loc[string]
    std = stats['std'].loc[string]
    plt_ax = plt.gca()
    plt_ax.text(.1, .6, "Mean: %.1f \n std: %.1f" % (mean,std), transform=plt_ax.transAxes)

def plot_mw(dataframe,foldername):
    "Plot the distribution of the molar weight"
    sns.set_palette(sns.color_palette("bright"))
    plt.subplots(figsize=(7.0,4.0))
    sns.histplot(dataframe.droplevel(0, axis=1), x='Mw', bins=50, stat='percent')
    plt.ylabel('Percent of dataset')
    plt.xlabel('Molar weight [Da]')
    plt.tight_layout()
    plt.savefig(foldername+'/'+'Mw'+'_dist.png',dpi=320)
    plt.close()

def plot_number_of_heavy_atoms(dataframe,foldername):
    "Plot the distribution of the heavy atoms"
    sns.color_palette("bright")
    plt.subplots(figsize=(7.0,4.5))
    sns.histplot(dataframe.droplevel(0, axis=1), x='# heavy atoms', bins=25, discrete=True)
    plt.ylabel('Number of molecules')
    plt.xlabel('Number of heavy atoms')
    plt.tight_layout()
    plt.savefig(foldername+'/'+'heavyAtoms'+'_dist.png',dpi=320)
    plt.close()

    tmp = dataframe[dataframe.columns[dataframe.columns.get_level_values(0)
    == 'Atoms'].intersection(dataframe.columns[dataframe.columns.get_level_values(1)
        != 'H'])].stack()
    tmp['type'] = tmp.index.get_level_values(1)
    element_plot = sns.catplot(tmp, x='type',y="Atoms", kind='bar')
    element_plot.set(xlabel="Element", ylabel="Count")
    plt.tight_layout()
    plt.savefig(foldername+'/'+'heavyAtoms'+'val_dist.png',dpi=320)
    plt.close()

def plot_atomic_ratios(dataframe, foldername):
    "Plot the atomic ratios"
    tmp = dataframe[dataframe.columns[dataframe.columns.get_level_values(0)
    == 'Atom ratios'].intersection(dataframe.columns[dataframe.columns.get_level_values(1)
        != 'H'])].stack()
    tmp['type'] = tmp.index.get_level_values(1)
    atom_ratio_plt = sns.catplot(tmp, x='type',y="Atom ratios", kind='box')
    atom_ratio_plt.set(xlabel="Ratio", ylabel="")
    atom_ratio_plt.set_titles(row_template="{row_name}")
    plt.tight_layout()
    plt.savefig(foldername+'/atomic_ratios.png')
    plt.close()
