"""script to reformat dataframes with string lists so that they have a compound column"""
import joblib as jl 
datasets={'Gecko':'/scratch/work/sandsth2/DataSets/gecko/gecko_full.dump',
         'Wang':'/scratch/work/sandsth2/DataSets/acp-17-7529-2017-supplement/wang_data.dump',
         'QM9':"/scratch/work/sandsth2/DataSets/QM9/qm9_data.dump",
         'Quinones':"/scratch/work/sandsth2/DataSets/quinones/Tabor/tabor_nosulf.dump",
         'ExpMoNA':"/scratch/work/sandsth2/DataSets/MassBank/expmona_cleaned.dump",
          'nablaDFT':"/scratch/work/sandsth2/DataSets/nablaDFT/nablaDFT.dump", 
         'IS':"/scratch/work/sandsth2/DataSets/Pesticides/IS_pesticides.dump",
         'RS':"/scratch/work/sandsth2/DataSets/Pesticides/RS_pesticides.dump"}
for key in datasets.keys():
    df = jl.load(datasets[key])
    df['compound'] = df.index
    df.to_csv(key+".csv")
