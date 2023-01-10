import pandas as pd
import os, sys
ROOT_DIR = os.path.dirname(os.path.abspath('/scratch/work/sandsth2/SCRIPTS/MolAnalysis'))
sys.path.append(ROOT_DIR)
from MolAnalysis.descriptors import CustomDescriptors
import joblib as jl 

MoNa_exp = pd.read_csv("MoNA-export-Experimental_Spectra.smiles.unique.-gt", header=None, on_bad_lines='skip', low_memory=False);
MoNa_exp.rename(columns = {0: "SMILES"}, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='(CC(=O)O3)O)C'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1) N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)O'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='N/ACCC1(C(=O)NCNC1=O)c2ccccc2'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='O=C(N(C(C=CC=C1)=C1C(N(C)[C@@]2([H])CC3=CC=CC=C3)=O)C2=N4)C5=C4C=CC=C6'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='O=C(N[C@@H](CCCCCC(CC)=O)C(N[C@@H](CC1=CN(OC)C2=C1C=CC=C2)C3=O)=O)[C@@H]4N(C([C@H]([C@H](CC)C)N3)=O)CCCC5'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='O=C([C@H](CC)C)O[C@H]1CCC=C2C1[C@@H](CC[C@@H](O)C[C@@H](O)CC(OC)=O)[C@@H](C)C=C3'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='O=C1C2=C(C=C(C)C=C2O)OC3=CC(O)=CC(C(OC)=O)=C32'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='C1=CC(C(OC)=O)=C(OC2=CC(C)=CC(O)=C2C(O)=O)C(OC)=C2'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='OC1=CC=C(CC(C(NC(C(CC)C)C(OC(C(CCCCCCCCCC)C)CC(NC(C(NC(C(NC(C(NC2CCC(N)=O)=O)C)=O)C)=O)C(O)C)=O)=O)=O)NC2=O)C=C2'].index, inplace=True)
MoNa_exp.drop(MoNa_exp.loc[MoNa_exp['SMILES']=='OC1=CC(C(OC)=O)=C(OC2=CC(C)=CC(O)=C2C(O)=O)C(OC)=C2'].index, inplace=True)
CustomDescriptors(MoNa_exp) 
jl.dump(MoNa_exp, "MoNa_exp.dump")
#MoNa_exp.columns