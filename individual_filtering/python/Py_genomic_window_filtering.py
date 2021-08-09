import pandas as pd
import sys

#script to calculate fraction of sites flagged as "bad" in sliding genomic windows (bad sites defined as those with allelic balance > 0.05 or <0.95). 

file_in = "/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/pileup_sep/"+sys.argv[1]+"_sep.txt"

print(file_in)

dat = pd.read_table(file_in, header=0,delim_whitespace=True)

#select chromosomes provided in command line argument
dat = dat.loc[(dat['V1'] == sys.argv[2] )

#select sex chromosomes scaffolds
#dat = dat.loc[(dat['V1'] == "Chr_Y_A" )  | (dat['V1'] == "Chr_Y_B" ) | (dat['V1'] == "scaffold_17" )| (dat['V1'] == "scaffold_18" )| (dat['V1'] == "scaffold_210" )| (dat['V1'] == "scaffold_227" )| (dat['V1'] == "scaffold_230" )| (dat['V1'] == "scaffold_240" )| (dat['V1'] == "scaffold_250" )| (dat['V1'] == "scaffold_277")| (dat['V1'] == "scaffold_497") ]     





#list unique set of scaffolds
scaffolds = pd.DataFrame(set(dat.V1))

#empty dataframe where results will be written
results = pd.DataFrame(columns=('chromosome','win_start','win_end', 'sumbad', 'window_length'))

for j in range(0,len(scaffolds)):
	df = dat.loc[dat['V1'] == scaffolds.iloc[j,0]]
	print(scaffolds.iloc[j,0])
	for i in range(0, df.iloc[(len(df)-1),1] - 99, 100):
		window= df.loc[(df['V2'] >= i) & (df['V2'] <= (i+350))]
		sumbad = window['bad'].sum()
		winlen = len(window.index)
		chrom = df.iloc[0,0]
		row = [chrom,i,(i+350),sumbad,winlen]
		results.loc[len(results)] = row
 
 
file_out = "/plas1/george.sandler/marchantia_popgen_newgenome/individual_filtering/windows_analysis_output/"+sys.argv[1]+"_"+sys.argv[2]+"_UV_windows.txt"
 
results.to_csv(file_out, header=None, index=None, sep='\t', mode='a')

exit()

