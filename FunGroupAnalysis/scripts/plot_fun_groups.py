import joblib as jl 
import seaborn as sns
import pandas as pd
import sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def main():
    fg_datasets={'IS':"simpol-outputs/IS_SIMPOLgroups_no-ring_nitrophenol.csv",
         'RS':"simpol-outputs/RS_SIMPOLgroups_no-ring_nitrophenol.csv"}
    threshold=0.00

    for key in fg_datasets:
        df = pd.read_csv(fg_datasets[key])
        df = process_data(df)
        df_c = count_values(df, threshold)
        plot_fun_group_hist(df,key)
        plot_fun_group_count_bar(df_c, threshold, key)  

if __name__ == "__main__":

    main()

def plot_fun_group_hist(df,key):
    """Plot histogram distribution of functional groups.
    df - pandas dataframe, key - name of dataset (str)"""
    plt.figure(figsize=(3,2),dpi=300)
    plt.rc('font', size=10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    df.sum(axis=1).hist(bins=20, ec='k')
    plt.ylabel('count')
    plt.xlabel('Number of functional groups')
    plt.grid(False)
    plt.subplots_adjust(bottom=0.4)
    plt.savefig(''.join(('images/'+key,'_','no-thres_FG_hist.png')), transparent=True, bbox_inches='tight')
    plt.close()
    
def process_data(df):
    """Remove unnecessary columns 1) zeroeth group?? 2) ester, all = ester,nitro+ester
    df - pandas dataframe"""
    df.set_index('compound', inplace=True)
    df.drop(['zeroeth group', 'ester, all'], axis=1, inplace=True)
    df.rename(columns={'carbon number on the acid-side of an amide (asa)':"#C on acid-side of amide"},inplace=True)
    return df[df.columns[df.sum() != 0]]

def count_values(df,threshold):
    """Count the different number of times each fun group is found in the molecules
    df - pandas dataframe, threshold - frequency threshold (0.0-1.0) for inclusion of functional group(float)"""
    lista=[]
    for column in df.columns: 
        lista.append(df[column].value_counts(normalize=True))
    df_c = pd.DataFrame(lista)
    df_c = pd.concat([df_c[[1,2,3]],df_c.drop([0,1,2,3],axis=1).sum(axis=1)],axis=1)
    df_c.columns = ['1','2','3',r'$\geq$4']
    df_c['total'] = df_c.sum(axis=1)
    return df_c[df_c['total'] >= threshold]

def plot_fun_group_count_bar(df_c, threshold, key):
    """Plot the count (normalized) and colored by different values as together with an insert of the fg histogram. df - pandas dataframe of value counts,threshold - frequency threshold (0.0-1.0) for inclusion of functional group(float), key - name of dataset (str)"""
    plt.figure(figsize=(7,9.5),dpi=300)
    plt.rc('font', size=10)
    ax = df_c.sort_values('total',ascending=False).drop('total',axis=1).plot(kind='barh', stacked=True)
    plt.legend(title='Occurence')
    fig = matplotlib.pyplot.gcf()
    plt.yticks(rotation=25)
    plt.subplots_adjust(left=0.3)
    plt.savefig(''.join(('images/'+key,'_', str(threshold), '_FG.png')),transparent=True, bbox_inches='tight')
    ax.text(0.4,10,'Activated fun. grps: '+str(len(df_c))+'\nThreshold = '+str(threshold))
    im = plt.imread(''.join(('images/'+key,'_','no-thres_FG_hist.png'))) # insert local path of the image.
    newax = fig.add_axes([0.35,0.5,0.37,0.37], anchor='NE', zorder=1)
    newax.imshow(im)
    newax.axis('off')
    plt.savefig(''.join(('images/'+key,'_insert_', str(threshold), '_FG.png')),transparent=True, bbox_inches='tight')
    plt.close()
    