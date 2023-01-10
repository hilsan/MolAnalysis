import seaborn as sns
import matplotlib.pyplot as plt
from  rdkit import Chem
from rdkit.Chem import (Descriptors, AllChem, rdchem, rdmolops)
import pandas as pd
import os
from collections import defaultdict
import numpy as np


class CustomDescriptors(pd.DataFrame):
    def __init__(self, *args, **kw):
        super(CustomDescriptors, self).__init__(*args, **kw)
        if  'mol' not in '\t'.join([str(x) for x in self.columns]) :
            self['mol'] = self['SMILES'].apply(lambda x: AllChem.MolFromSmiles(x)  if(np.all(pd.notnull(x))) else np.NaN)
        if 'Mw' not in '\t'.join([str(x) for x in self.columns]):
                self['Mw'] = self['mol'].apply(lambda x: Descriptors.MolWt(x) if(np.all(pd.notnull(x))) else np.NaN)
        if 'Func. groups' not in self.columns: 
            #Attach the functional groups
            tmp = self['mol'].apply(lambda x: countFunctionalGroups(x) if(np.all(pd.notnull(x))) else np.NaN).apply(pd.Series)
            tmp.columns = pd.MultiIndex.from_product([['Func. groups'],tmp.columns])
            self[tmp.columns] = pd.DataFrame(tmp)
        if 'Atom ratios' not in self.columns:
            #Attach the atomic composition and ratios
            tmp = self['mol'].apply(lambda x: getAtomicComposition(x) if(np.all(pd.notnull(x))) else np.NaN)
            self['comp_hash'] = tmp
            tmp2 = tmp.apply(pd.Series)
            self['# heavy atoms'] = tmp2[tmp2.columns[tmp2.columns.str.contains('H') == False]].sum(axis=1)
            tmp2.columns = pd.MultiIndex.from_product([['Atoms'],tmp2.columns])
            self[tmp2.columns] = pd.DataFrame(tmp2)
            tmp = self['comp_hash'].apply(lambda x: getAtomicRatios(x) if(np.all(pd.notnull(x))) else np.NaN).apply(pd.Series)
            tmp.columns = pd.MultiIndex.from_product([['Atom ratios'],tmp.columns])
            self[tmp.columns] = pd.DataFrame(tmp)
        li = []
        for item in self.columns:
            if type(item) == str: 
                li.append(("",item))
            else: 
                li.append(item)
        self.columns = pd.MultiIndex.from_tuples(li)


def countFunctionalGroups(mol):
    fg_list = ["C=O", "C#N", "CN", "C#C"]
    res = {}
    for group in fg_list: 
        functional_group = Chem.MolFromSmiles(group)
        matches = mol.GetSubstructMatches(functional_group)
        res[group] = len(matches)
    res['# of func. groups'] = sum(res.values())
    return res

def getMolFromSmiles(smiles):
    if isinstance(smiles, str):
        mol = Chem.MolFromSmiles(smiles)
    else: 
        mol = np.NaN
        print('smiles'+smiles+'is not valid')
    return mol

def getAtomicComposition(mol):
    if mol:
        # Add hydrogen atoms--RDKit excludes them by default
        molecule_with_Hs = Chem.AddHs(mol)
    comp = defaultdict(int)
    # Get atom counts
    for atom in molecule_with_Hs.GetAtoms():
        comp[atom.GetSymbol()] += 1
    #comp['# heavy atoms'] = sum(comp.values()) - comp['H']
    
    # If charged, add charge as "atomic number" 0
    # charge = rdmolops.GetFormalCharge(molecule_with_Hs)
    #comp['charge'] = charge
    return comp

def getAtomicRatios(comp):
    ratios = {}
    if comp['C']:
        ratios['O:C'] = comp['O']/comp['C']
        ratios['H:C'] = comp['H']/comp['C']
        ratios['N:C'] = comp['N']/comp['C']
    else: 
        ratios['O:C'] = np.NaN
        ratios['H:C'] = np.NaN
        ratios['N:C'] = np.NaN
    return ratios