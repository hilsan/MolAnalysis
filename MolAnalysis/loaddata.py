
import joblib as jl 
import pandas as pd


def loadQM9(): 
    path = "/scratch/work/sandsth2/DataSets/" +"QM9/qm9_data.dump"
    return jl.load(path) 

def loadGecko(): 
    path = "/scratch/work/sandsth2/DataSets/" +"gecko/gecko_dataset.dump"
    return jl.load(path) 

def loadWang():
    path = "/scratch/work/sandsth2/DataSets/"+"acp-17-7529-2017-supplement/wang_data.dump"
    return jl.load(path)

def loadPesticides(): 
    RS_name = "/scratch/work/sandsth2/DataSets/Pesticides/RSdata.csv"
    IS_name = "/scratch/work/sandsth2/DataSets/Pesticides/ISdata.csv"
    return pd.read_csv(IS_name), pd.read_csv(RS_name)

def loadNablaDFT(): 
    path = "/scratch/work/sandsth2/DataSets/nablaDFT/nablaDFT.dump"
    return pd.read_csv(path)

def loadQuinones(): 
    path = "/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_nosulf.dump"
    #tabor_all: "/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_all.dump"
    return jl.load(path)

def loadExpMoNA(): 
    path = "/Users/sandsth2/triton-workdir/DataSets/MassBank/MoNa_exp.dump"
    return jl.load(path)

def loadInSilMoNA(): 
    path ="/Users/sandsth2/triton-workdir/DataSets/MassBank/MoNa_silico.dump"
    return jl.load(path)