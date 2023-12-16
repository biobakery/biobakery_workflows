import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--table", help="phylophlan table", type=str, required=True)
parser.add_argument("--output", help="output file", type=str, required=True)
args = parser.parse_args()

phylophlan = pd.read_csv(args.table, skiprows=[0,1,2], sep='\t')
phylophlan.columns = ['mag', 'sgb', 'ggb', 'fgb', 'ref']
phylophlan['taxon'] = ""
phylophlan['dist'] = 1

phylophlan = phylophlan.reset_index()
for i, row in phylophlan.iterrows():
	if float(row['sgb'].split(":")[3]) < 0.05:
		phylophlan.iloc[i, phylophlan.columns.get_loc('taxon')] = row['sgb'].split(":")[2]
		phylophlan.iloc[i, phylophlan.columns.get_loc('dist')] = float(row['sgb'].split(":")[3])
	elif float(row['ggb'].split(":")[3]) < 0.15:
		phylophlan.iloc[i, phylophlan.columns.get_loc('taxon')] = row['ggb'].split(":")[2]
		phylophlan.iloc[i, phylophlan.columns.get_loc('dist')] = float(row['ggb'].split(":")[3])
	elif float(row['fgb'].split(":")[3]) < 0.3:
		phylophlan.iloc[i, phylophlan.columns.get_loc('taxon')] = row['fgb'].split(":")[2]
		phylophlan.iloc[i, phylophlan.columns.get_loc('dist')] = float(row['fgb'].split(":")[3])
	else:
		phylophlan.iloc[i, phylophlan.columns.get_loc('taxon')] = "UNKNOWN"
		phylophlan.iloc[i, phylophlan.columns.get_loc('dist')] = 1

phylophlan.to_csv(args.output, sep='\t', index=False)
