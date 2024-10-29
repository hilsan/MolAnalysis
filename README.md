# A package to gather molecular representation and structural properties

Version of this project used for atmospheric compound similarity analysis published here:  DOI:10.5281/zenodo.14007835. 
This project contains functions to create molecular fingerprints, atomic ratios and molecular size from RDKit for datasets of SMILES strings. 

# SCRIPTS:
 descriptors.py - contains a class with functions to calculate atomic ratios, functional group statistics, elemental composition, molar weight, number of heavy atoms and atomic ratios for a dataframe with a column named 'SMILES
 compile_analysis.py - function to organize descriptor statistics and plot them.
 fingerprints.py - Fingerprints: a class to make molecular fingerprints for a whole list of SMILES strings in a column of a dataframe called 'SMILES'.  MACCSAnalysis: a class to gather statistics and correlations between features in MACCS fingerprints in dataframe.
