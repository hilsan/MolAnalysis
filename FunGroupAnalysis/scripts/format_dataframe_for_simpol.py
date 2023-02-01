"""script to reformat dataframes with string lists so that they have a compound column"""
import joblib as jl 
datasets={'IS':"/scratch/work/sandsth2/DataSets/Pesticides/IS_pesticides.dump",
         'RS':"/scratch/work/sandsth2/DataSets/Pesticides/RS_pesticides.dump"}
for key in datasets.keys():
    df = jl.load(datasets[key])
    df['compound'] = df.index
    df.to_csv(key+".csv")
