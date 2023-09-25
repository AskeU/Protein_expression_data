# -*- coding: utf-8 -*-
"""
Created on Mon May  2 09:11:27 2022

@author: askung
"""

####Import packages####
import pandas as pd
import os
from Bio import SeqIO
import glob
from Bio.Blast.Applications import NcbiblastnCommandline
'''
directory_path = "C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/NGS/2023_08_16_Calm_mem_2/"  # Replace with your directory path
database="C:/blast/BLAST_DBs/calmodulin/calm_db"
output_path=f'{directory_path}/Calmodulin' #Where should the results go
fastqpath = f'{directory_path}/pass'
'''

################BLASTER##############

#Main script calling above
def blaster(database,output_path,fastqpath):
    
    fastq_files = [os.path.join(fastqpath, file) for file in os.listdir(fastqpath) if file.endswith(".fastq")]
    output_fasta = "output.fasta"
    print("fastqfiles_found")
    #Make fastQ into fasta
    with open(output_fasta, "w") as output_handle:
        for fastq_file in fastq_files:
            sequences = SeqIO.parse(fastq_file, "fastq")
            SeqIO.write(sequences, output_handle, "fasta")
    print("fasta_files created")
    blastn = NcbiblastnCommandline(cmd="C:/blast/bin/blastn.exe", query=output_fasta, db=database, outfmt=6)
    stdout, stderr = blastn()
    blast_results=stdout
    print("blast_done")
    # Convert BLAST results to CSV format
    csv_header = "read_id,alignment_genome,pident,length,mismatch,gapopen,alignment_strand_start,alignment_strand_end,alignment_genome_start,alignment_genome_end,evalue,bitscore"
    csv_content = blast_results.replace('\t', ',')
    
    # Save BLAST results to a file
    with open(f"{output_path}/blast_results.csv", "w") as result_file:
        result_file.write(csv_header + '\n' + csv_content)

    print("BLAST results saved to 'blast_results.csv'")

#blaster(database,output_path,fastqpath)
##############Merger####################

#alignment_data=output_path

def Merger(Alignment_Data,FastQ_path): 
    os.chdir(Alignment_Data)
    Summary = pd.read_csv('blast_results.csv',sep=",")
    
    os.chdir(FastQ_path)
    Filelist =(glob.glob("*.fastq"))
    
    SeqData = []
    for i in Filelist:
        for record in SeqIO.parse(i, "fastq"):
            SeqData.append([record.id,record.seq])
    
    SeqData = pd.DataFrame(data=SeqData)
    
    SeqData.columns=["read_id","Sequence"]
    
    Merged_Data= pd.merge(SeqData,Summary,on=["read_id"])
    
    Merged_Data["alignment_direction"]= "+"
    Merged_Data["alignment_direction"][Merged_Data["alignment_genome_start"]>Merged_Data["alignment_genome_end"]]="-"
    
    #Swap start, end at the points where they are reversed
    condition = Merged_Data["alignment_direction"] == "-"
    Merged_Data.loc[condition, ['alignment_genome_start', 'alignment_genome_end']] = Merged_Data.loc[condition, ['alignment_genome_end', 'alignment_genome_start']].values
    return Merged_Data



#merged_data = Merger(alignment_data, fastqpath)

