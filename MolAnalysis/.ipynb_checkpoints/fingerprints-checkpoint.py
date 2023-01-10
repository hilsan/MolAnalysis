from rdkit import (Chem, DataStructs)
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
import sys
import joblib as jl 
import opendatasets as od
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
import seaborn as sns

class fingerprints: 
    #Default values set from the Lumiaro Paper 
    def __init__(self, 
                fpType = None, 
                r = 2, 
                n_bits=2048, 
                n_bits_per_hash = 16, 
                max_path = 8,
                fp_size = 8192):

                self.fpType = fpType
                self.radius = r
                self.nBits = n_bits
                self.nBitsPerHash = n_bits_per_hash
                self.maxPath = max_path
                self.fpSize = fp_size

    def getFingerPrint(self,smiles, **kwargs):
        try:
            mol = AllChem.MolFromSmiles(smiles)
            if self.fpType == 'Morgan':
                fingerprint =  AllChem.GetMorganFingerprintAsBitVect(mol, radius = self.radius, nBits = self.nBits, **kwargs) 
            elif self.fpType == 'MACCS':
                fingerprint = MACCSkeys.GenMACCSKeys(mol, **kwargs)
            elif self.fpType == 'topological':    
                fingerprint = Chem.RDKFingerprint(mol, maxPath = self.maxPath, 
                                                  nBitsPerHash = self.nBitsPerHash, fpSize = self.fpSize, **kwargs)
            elif self.fpType == None:
                raise Exception('Error: Enter fingerprint type')
        except:
            fingerprint=np.NaN
            
        return fingerprint  


class MACCSAnalysis(pd.DataFrame):
    
        def __init__(self, *args, **kw):
            super(MACCSAnalysis, self).__init__(*args, **kw)
            self['MACCS'] = self['SMILES'].apply(lambda x: fingerprints(fpType='MACCS').getFingerPrint(x) if(np.all(pd.notnull(x))) else np.NaN)
            sns.set_palette(sns.color_palette("bright"))
            
        def maccsToArray(self,fp): 
            try:
                 bitfp =  np.frombuffer(fp.ToBitString().encode(), 'u1') - ord('0')
            except: 
                bitfp = np.NaN
            return bitfp
        
        def makeMACCSDF(self):
            return pd.DataFrame(self['MACCS'].apply(lambda x: self.maccsToArray(x) if(np.all(pd.notnull(x))) else np.NaN).to_list())
        
        def get_counts(self):
            feature_counts = self.makeMACCSDF().sum()
            nz_features = feature_counts[feature_counts != 0]/len(self.makeMACCSDF())
            return feature_counts, nz_features
        
        def shared_features(self, nz_features, percentages):
            n_shared_features=[]
            for i in percentages:
                n_shared_features.append(len(nz_features[nz_features > i]))
            return n_shared_features
        
        def group_correlated_features(self, nz_features, threshold, count_threshold): 
            nz_corr = self.makeMACCSDF()[nz_features.index.to_list()].corr()
            a = nz_corr.where((nz_corr  > threshold) | (nz_corr  < -threshold), np.NaN)
            corr_features = []
            for i in a.columns: 
               corr_features.append(list(np.where(a.loc[i] > threshold))[0])
            return pd.Series(corr_features).apply(lambda x: len(x))
            
        def plot_heatmap(self, nz_features, filename, threshold):
            nz_corr = self.makeMACCSDF()[nz_features.index.to_list()].corr()
            a = nz_corr.where((nz_corr  > threshold) | (nz_corr  < -threshold), np.NaN)
            heatplot = sns.heatmap(a, cmap="Blues", xticklabels=nz_corr.columns, yticklabels=nz_corr.columns)
            fig = heatplot.get_figure()
            fig.savefig(filename+"_heatmap.png") 

        def plot_feature_distribution(self, shared_features, nz_features, percentages, filename): 
            #Plotta hur många  feature finns i mer än x procent av datasetet)
            fig, axs = plt.subplots(1,2,dpi=200)
            fig.subplots_adjust(wspace=0.2, hspace=0.05)
            axs[0].plot(percentages, shared_features, marker='o', markersize=2.5, linewidth=1)
            axs[0].set_ylabel('Number of features')
            axs[1].hist(nz_features, 20,facecolor='b', histtype='bar',rwidth=3.0)
            plt.setp(axs, xlim=[-0.04,1], xlabel = 'Percent of dataset');
            plt.savefig(filename+"_fdist.png")

        def plot_FeatureFP(self, feature_counts, filename): 
            fig, axs  = plt.subplots(figsize=[7,3.6],dpi=200)
            axs.bar(np.arange(167),  feature_counts/len(self))
            plt.setp(axs, xlim=[-0.04,166], ylim=[0,1], xlabel = 'Feature ID', ylabel='Feature count');
            plt.savefig(filename+"_featurefp.png")
