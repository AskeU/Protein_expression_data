# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 13:51:55 2022

@author: askung
"""

from Bio.Seq import Seq
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
import matplotlib.colors
import math
import glob
import os
import re
import importlib.util
import matplotlib.colors as mcolors
path = r'C:\Users\askung\OneDrive - Danmarks Tekniske Universitet\WP1\Article writing\Article_Course\2022_12_06_NGS_Hamilton_Script_package'
os.chdir(path)
print(os.getcwd())

###Import scripts
from Double_barcode_Demultiplexer_vectorized import NGS_Demultiplexer_vectorized
from Demultiplexed_To_Stuctured_With_Lum_data import demultiplex_to_lum
from Merger_blast import blaster
from Merger_blast import Merger

os.chdir(r"C:\Users\askung")

####General data####
#Import primers (For barcode positions)
Fw_Primers = pd.read_excel ('C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/Barcode_Analysis/Fw_Primers.xlsx')
Rev_Primers = pd.read_excel ('C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/Barcode_Analysis/Rev_Primers.xlsx')

#################
######Paths######
#Blast: 
database="C:/blast/BLAST_DBs/tf3/tf3_db"


#NGS_data general path
NGS_path = "C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/"

#Fill in "First_primers_used (For NGS), ORF_length and Alignment_Data
#Path to find merger data 

Alignment_Data=f'{NGS_path}/2023_07_26_spa_tf3/TF3_alignment'
Path_With_CSVname=f'{Alignment_Data}Merged_Nanopore.csv'

#Which was the first primer for NGS?
First_primer_used = 25 #(For NGS data)
#How long was the sequence?
ORF_Length = 885




#Path to find luminescence data from Hamilton (Different path than the rest)
luminescence_path = "C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/Automation_Work/2023_07_07_SpA_TF3"

#FastQ files (For merger)--> NOT BAM files
FastQ_path = f'{NGS_path}/2023_07_26_spa_tf3'


#########Script-paths below, not needed to fill########


#Merger combines the data and finds the alignment direction
Merger_path = f'{path}Nanopore_Merger.py'
#Path for running "Demultiplex_script"
Demultiplex_Path = f'{path}Double_barcode_Demultiplexer.py'
#Path for running "After_Demultiplex_script"
AfterDemultiplex_Path = f'{path}Demultiplexed_To_Stuctured_With_Lum_data.py'
#Outputs for data analysis

#main_blaster blasts the fastQ´s vs the "genome". Create them by running 
#C:/blast/bin/makeblastdb.exe -in "C:/blast/BLAST_DBs/[fasta_genome].fa" -dbtype nucl -out [genome]_db
blaster(database,Alignment_Data,FastQ_path)

#Run the demultiplexer (Not needed if only Hamilton)
merged = Merger(Alignment_Data,FastQ_path)
merged.to_csv(Path_With_CSVname, index = False)
print("Done with merger")
#Run the demultiplexer


#Alignment_merged (for Demultiplex)
Sequence_samples = pd.read_csv(Path_With_CSVname)
print('\n### Shape of seq samples:', Sequence_samples.shape)
after_demultiplex = NGS_Demultiplexer_vectorized(Fw_Primers,Rev_Primers,Sequence_samples)

#Remove contaminations, if this is removed, uncertain results will not be filtered
after_demultiplex_qc = after_demultiplex[after_demultiplex["decontamination"]==True]

after_demultiplex_qc.to_csv(f'{Alignment_Data}/after_demultiplex.csv', index = False)
print("Done with demultiplex and QC")

#Run the secondary data-crunching with luminescence-data

NGS_Filtered = demultiplex_to_lum(Alignment_Data, luminescence_path, First_primer_used, ORF_Length)
print("Done with script")
