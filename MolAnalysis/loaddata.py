
"Contains functions which load different datasets"

import joblib as jl
import pandas as pd


def load_qm9():
    "Loads a qm9 dataset"
    path = "/scratch/work/sandsth2/DataSets/" +"QM9/qm9_data.dump"
    return jl.load(path)

def load_gecko():
    "Loads a gecko dataset"
    path = "/scratch/work/sandsth2/DataSets/" +"gecko/gecko_full.dump"
    return jl.load(path)

def load_wang():
    "Loads a wang dataset"
    path = "/scratch/work/sandsth2/DataSets/"+"acp-17-7529-2017-supplement/wang_data.dump"
    return jl.load(path)

def load_pesticides():
    "loads a pesticides dataset"
    rs_name = "/scratch/work/sandsth2/DataSets/Pesticides/RSdata.csv"
    is_name = "/scratch/work/sandsth2/DataSets/Pesticides/ISdata.csv"
    return pd.read_csv(is_name), pd.read_csv(rs_name)

def load_nabladft():
    "loads a nabladft dataset"
    path = "/scratch/work/sandsth2/DataSets/nablaDFT/nablaDFT.dump"
    return pd.read_csv(path)

def load_quinones():
    "loads a quinones dataset"
    path = "/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_nosulf.dump"
    #tabor_all: "/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_all.dump"
    return jl.load(path)

def load_exp_mona():
    "loads a exp. mona dataset"
    path = "/Users/sandsth2/triton-workdir/DataSets/MassBank/MoNa_exp.dump"
    return jl.load(path)

def load_insil_mona():
    "loads a in silico mona dataset"
    path ="/Users/sandsth2/triton-workdir/DataSets/MassBank/MoNa_silico.dump"
    return jl.load(path)
