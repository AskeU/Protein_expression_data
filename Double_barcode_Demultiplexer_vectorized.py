# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:06:52 2022

This script runs a small 9-barcode twice for each

@author: askung
"""

from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors



def NGS_Demultiplexer_vectorized(Fw,Rev,Sequences): 
    #Put in plate and position number 
    print('# Shape of Sequence_samples:', Sequences.shape)

    Sequences["PlateNo"]=0
    Sequences["PosNo"]=0
    #List to save all the barcodes found
    Sequences["full_fw_list"]='A'
    Sequences["full_rev_list"]='A'
    

    Sequence_samplesMinus = Sequences[Sequences["alignment_direction"].str.contains("-")].reset_index()
    Sequence_samplesPlus = Sequences[-(Sequences["alignment_direction"].str.contains("-"))].reset_index()
    

    #Find barcode possibilities:
    Fw_Barcodes = Fw.iloc[:,2]
    Rev_Barcodes = Rev.iloc[:,2]
    
    #Barcode runs for Fw/rev and reverse complement of them

    #For all of the Fw barcodes--> Find the ones that contains the barcodes
    for i in range(Fw_Barcodes.shape[0]):
        #Defines the 9 base barcode
        Barcode = Seq(Fw_Barcodes[i][0:9])
    
        Barcode = str(Barcode)
        #Positions samples based on tags:
        for j in range(len(Sequence_samplesPlus["Sequence"])):
            First_Barcode = str(Sequence_samplesPlus["Sequence"][j]).find(Barcode)
            #First Barcode is found between 0 and the start position of the alignment
            if First_Barcode>=0 and First_Barcode<=Sequence_samplesPlus["alignment_strand_start"][j]:
                Second_Barcode = (str(Sequence_samplesPlus["Sequence"][j][First_Barcode+1:First_Barcode+16]).find(Barcode))
                if Second_Barcode>=0:
                    Sequence_samplesPlus["PosNo"][j]= i+1
                    Sequence_samplesPlus["full_fw_list"][j] = f'{Sequence_samplesPlus["full_fw_list"][j],i+1}'
  
    print("Done with barcode_fw")

    for i in range(Fw_Barcodes.shape[0]):
        #Defines the 9 base barcode
        Barcode = Seq(Fw_Barcodes[i][0:9]).reverse_complement()
        Barcode = str(Barcode)
        #Positions samples based on tags:
        for j in range(len(Sequence_samplesMinus["Sequence"])):
            First_Barcode = str(Sequence_samplesMinus["Sequence"][j]).find(Barcode)
    #If the first barcode is above the alignment start (It is reverse complement so opposite here)
            if First_Barcode>=Sequence_samplesMinus["alignment_strand_start"][j]:
                Second_Barcode = (str(Sequence_samplesMinus["Sequence"][j][First_Barcode+1:First_Barcode+16]).find(Barcode))
                if Second_Barcode>=0:
                    Sequence_samplesMinus["PosNo"][j]= i+1
                    Sequence_samplesMinus["full_fw_list"][j] = f'{Sequence_samplesMinus["full_fw_list"][j],i+1}'

    print("Done with barcode_fw 2")
    
                    
    for i in range(Rev_Barcodes.shape[0]):
        #Defines the 9 base barcode
        Barcode = Seq(Rev_Barcodes[i][0:9]).reverse_complement()
        Barcode = str(Barcode)
        #Positions samples based on tags:
        for j in range(len(Sequence_samplesPlus["Sequence"])):
            First_Barcode = str(Sequence_samplesPlus["Sequence"][j]).find(Barcode)
            #If the barcode is above the end of the alignment:
            if First_Barcode>=Sequence_samplesPlus["alignment_strand_end"][j]:
                Second_Barcode = (str(Sequence_samplesPlus["Sequence"][j][First_Barcode+1:First_Barcode+16]).find(Barcode))
                if Second_Barcode>=0:
                    Sequence_samplesPlus["PlateNo"][j]= i+1
                    Sequence_samplesPlus["full_rev_list"][j] = f'{Sequence_samplesPlus["full_rev_list"][j],i+1}'

    print("Done with barcode_rev")
    
    for i in range(Rev_Barcodes.shape[0]):
        #Defines the 9 base barcode
        Barcode = Seq(Rev_Barcodes[i][0:9])
        Barcode = str(Barcode)
        #Positions samples based on tags:
        for j in range(len(Sequence_samplesMinus["Sequence"])):
            First_Barcode = str(Sequence_samplesMinus["Sequence"][j]).find(Barcode)
    #If the barcode is before the end position (which is the full seq - end position)
            if First_Barcode<=len(Sequence_samplesPlus["Sequence"])-Sequence_samplesPlus["alignment_strand_end"][j]:
                Second_Barcode = (str(Sequence_samplesMinus["Sequence"][j][First_Barcode+1:First_Barcode+16]).find(Barcode))
                if Second_Barcode>=0:
                    Sequence_samplesMinus["PlateNo"][j]= i+1
                    Sequence_samplesMinus["full_rev_list"][j] = f'{Sequence_samplesMinus["full_rev_list"][j],i+1}'
    
    print("Done with barcode_rev2")
                   
    Sequence_samples_Demultiplexed=Sequence_samplesMinus.append(Sequence_samplesPlus)
    

    Test_Info= Sequence_samples_Demultiplexed[Sequence_samples_Demultiplexed["PosNo"]>0]
    Test_Info2=Test_Info[Test_Info["PlateNo"]>0].reset_index()
    
    Test_Info2[Test_Info2["alignment_genome"]!="*"] #Remove non-aligned pieces
    
    Test_Info2["Positions"] = Test_Info2["PlateNo"].astype(str)+","+Test_Info2["PosNo"].astype(str)
    Test_Info2["Duplicates"]=0
    
    print("Starting finding duplicates")

    #Find the amount of duplicates in the set
    for i in range(len(Test_Info2["Duplicates"])):
        Test_Info2["Duplicates"][Test_Info2["Positions"]==Test_Info2["Positions"][i]] +=1
        
    
    #Change Test_Info2 to only include sequences with +3 duplicates:
    
    #Remove unnecessary stuff from lists with all the predicted barcodes
    Test_Info2["full_fw_list"]=Test_Info2["full_fw_list"].str.replace(r"(","").str.replace(r")","").str.replace(r"'","").str.replace(r"A, ","").str.replace(r'"',"")
    Test_Info2["full_rev_list"]=Test_Info2["full_rev_list"].str.replace(r"(","").str.replace(r")","").str.replace(r"'","").str.replace(r"A, ","").str.replace(r'"',"")

    #Find right positions of the ones with multiple hits
    Test_Info2= Test_Info2[Test_Info2["full_fw_list"].str.find("A")<0]
    Test_Info2= Test_Info2[Test_Info2["full_rev_list"].str.find("A")<0]

    Test_Info2=Test_Info2[Test_Info2["Duplicates"]>5]

    #Dataframe without the multihits:
    Test_Info3 = Test_Info2[Test_Info2["full_fw_list"].str.find(",")<0][Test_Info2["full_rev_list"].str.find(",")<0]
    
    #Extract the positions and find the ones that fit with the same start/end position
    #Data to find out what is going on:
    
    Test_Info_8back = Test_Info3.copy()[["full_fw_list", "full_rev_list","alignment_genome_start","alignment_genome_end","Positions","read_id"]].reset_index()
    Test_Info_8back["full_fw_list"] = Test_Info_8back["full_fw_list"].astype(int)-8
    #Go through the fw. primers 8 back to see if there has been contaminations
    print("Contamination remover")
    # df1 is the positions 8 behind the specific position, df2 is the original positions
    df1 = Test_Info_8back
    df2 = Test_Info3.copy()[["full_fw_list", "full_rev_list","alignment_genome_start","alignment_genome_end","Positions"]].drop_duplicates().reset_index()


    print("length with duplicates: ",len(df1))
    
    df2.full_fw_list = df2.full_fw_list.astype(int)
    lenDF2=len(df2)
    print("length without duplicates: ",lenDF2)


    df1['Index'] = df1.index
    df2['Index'] = df2.index
    
    #It now runs a decontamination s

    # Total number of rows in df2 (used for following progress)
    total_rows = len(df2)
    
    # Empty list to hold matched indices from df1
    match_idx1 = []
    
    # Go through the original file and compare it to the positions 8 back. 
    #If they are the same, there must be a contamination and it will be removed. 
    for idx, (idx2, row2) in enumerate(df2.iterrows()):
        
        # Print progress
        print(f"Processing row {idx + 1} of {total_rows}...")
        
        # Filter df1 to only include potential matches for the current row in df2
        potential_matches = df1[
            (df1['alignment_genome_start'] - 6 <= row2.alignment_genome_start) & 
            (df1['alignment_genome_start'] + 6 >= row2.alignment_genome_start)
        ]
        
        # Create a 'key' column to facilitate merging
        potential_matches['key'] = 1
        df2['key'] = 1
    
        # Merge potential_matches with the current row in df2
        merged = pd.merge(potential_matches, pd.DataFrame([row2]))
    
        # Filter the merged dataframe based on the conditions
        mask = (
            (merged['full_fw_list'] == merged['full_fw_list']) & 
            (merged['full_rev_list'] == merged['full_rev_list']) &
            (abs(merged['alignment_genome_start'] - merged['alignment_genome_start']) <= 6) & 
            (abs(merged['alignment_genome_end'] - merged['alignment_genome_end']) <= 6)
        )
    
        matched_indices = merged[mask]['Index'].unique()
        match_idx1.extend(matched_indices)


    # See how many contaminations were found: 
    print('Matches:', len(match_idx1), "len(match_idx2), Not needed")
    #Now we have an index for the samples 
    print("QC_File maker:")
    #Run through Test_Info3 and find the ones that should be removed
    Test_Info3["decontamination"] = True
    for i in df1.loc[match_idx1]["read_id"]:
        Test_Info3["decontamination"][Test_Info3["read_id"]==i] = False

    print("########Demultiplexing Done##########")
    return Test_Info3
 