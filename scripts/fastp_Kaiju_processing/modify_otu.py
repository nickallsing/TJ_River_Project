###Author: Nicholas Allsing
###Purpose: Prepare OTU Table for Analysis in R
###Modify OTU table using pandas
#The OTU table now needs to be modified in order to effectively analyze the data
import pandas
import os
#read in otu table
metagenome = pandas.read_csv('../../data/fastp_Kaiju_processing/output/kaiju_out/OTU_metagenome.csv')
#set as dataframe
dfmg = pandas.DataFrame(metagenome)
#drop first 7 columns
dfmg.drop(dfmg.columns[[0,1,2,3,4,5,6]], axis=1, inplace=True)
#rename species to "ID" and add last 3 row names
dfmg.iloc[-1, dfmg.columns.get_loc('species')] = "unclassified"
dfmg.iloc[-2, dfmg.columns.get_loc('species')] = "cannot be assigned to a (non-viral) species"
dfmg.iloc[-3, dfmg.columns.get_loc('species')] = "Viruses"
dfmg = dfmg.rename(columns={"species" : "ID"})
#write out new metagenome file
dfmg.to_csv('../../data/fastp_Kaiju_processing/output/metagenomemodified.csv', index=False)
