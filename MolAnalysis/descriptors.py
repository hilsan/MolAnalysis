"Module to add molecular descriptors to a dataframe containing a column of SMILES strings"
from collections import defaultdict
import pandas as pd
import numpy as np
from rdkit.Chem import (Descriptors, MolFromSmiles, AddHs)

class CustomDescriptors(pd.DataFrame):
    "Class with custom molecular descriptors to a dataframe containing a column of SMILES strings"
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        if  'mol' not in '\t'.join([str(x) for x in self.columns]) :
            self['mol'] = self['SMILES'].apply(lambda x: MolFromSmiles(x)
                if(np.all(pd.notnull(x))) else np.NaN)
        if 'Mw' not in '\t'.join([str(x) for x in self.columns]):
            self['Mw'] = self['mol'].apply(lambda x: Descriptors.MolWt(x)
                if(np.all(pd.notnull(x))) else np.NaN)
        if 'Func. groups' not in self.columns:
            #Attach the functional groups
            tmp = self['mol'].apply(lambda x: count_functional_groups(x)
                if(np.all(pd.notnull(x))) else np.NaN).apply(pd.Series)
            tmp.columns = pd.MultiIndex.from_product([['Func. groups'],tmp.columns])
            self[tmp.columns] = pd.DataFrame(tmp)
        if 'Atom ratios' not in self.columns:
            #Attach the atomic composition and ratios
            tmp = self['mol'].apply(lambda x: get_atomic_composition(x)
                if(np.all(pd.notnull(x))) else np.NaN)
            self['comp_hash'] = tmp
            tmp2 = tmp.apply(pd.Series)
            self['# heavy atoms'] = tmp2[tmp2.columns[tmp2.columns.str.contains('H')
            is False]].sum(axis=1)
            tmp2.columns = pd.MultiIndex.from_product([['Atoms'],tmp2.columns])
            self[tmp2.columns] = pd.DataFrame(tmp2)
            tmp = self['comp_hash'].apply(lambda x: get_atomic_ratios(x)
                if(np.all(pd.notnull(x))) else np.NaN).apply(pd.Series)
            tmp.columns = pd.MultiIndex.from_product([['Atom ratios'],tmp.columns])
            self[tmp.columns] = pd.DataFrame(tmp)
        tmp_list = []
        for item in self.columns:
            if isinstance(item, str):
                tmp_list.append(("",item))
            else:
                tmp_list.append(item)
        self.columns = pd.MultiIndex.from_tuples(tmp_list)


def count_functional_groups(mol):
    "Function to count the functional groups of a molecule. Very basic SMILES matching"
    fg_list = ["C=O", "C#N", "CN", "C#C"]
    res = {}
    for group in fg_list:
        functional_group = MolFromSmiles(group)
        matches = mol.GetSubstructMatches(functional_group)
        res[group] = len(matches)
    res['# of func. groups'] = sum(res.values())
    return res

def get_mol_from_smiles(smiles):
    "Function to convert a SMILES string to a molecule"
    if isinstance(smiles, str):
        mol = MolFromSmiles(smiles)
    else:
        mol = np.NaN
        print('smiles'+smiles+'is not valid')
    return mol

def get_atomic_composition(mol):
    "Function to get all element types in a mol struct"
    if mol:
        # Add hydrogen atoms--RDKit excludes them by default
        molecule_with_hs = AddHs(mol)
    comp = defaultdict(int)
    # Get atom counts
    for atom in molecule_with_hs.GetAtoms():
        comp[atom.GetSymbol()] += 1
    #comp['# heavy atoms'] = sum(comp.values()) - comp['H']
    # If charged, add charge as "atomic number" 0
    # charge = rdmolops.GetFormalCharge(molecule_with_Hs)
    #comp['charge'] = charge
    return comp

def get_atomic_ratios(comp):
    "Function to get atomic ratios"
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
