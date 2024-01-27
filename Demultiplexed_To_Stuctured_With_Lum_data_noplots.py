# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 09:14:54 2022

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
import matplotlib.colors as mcolors
import matplotlib.cm as cm


def demultiplex_to_lum(Output_Path, luminescence_path, First_primer_used,ORF_Length):
    #Read data from demultiplexer
    After = pd.read_csv(f'{Output_Path}/after_demultiplex.csv',sep=",")
    
    #Fix up data
    Seq_Data=After.sort_values("Positions")
    #Remove non-aligned
    Seq_Data= Seq_Data[Seq_Data["alignment_genome"]!="*"]

    #Find all plate positions
    Positions = [Seq_Data.Positions.unique()]
    
    #Array to fill data about probable position into, added a 0-array to start
    Combined_NGS = np.array([0,0,0,0,0,0,0,0])

    #Loop that goes through all the found plate-positions 
    for platepositions in range(len(Positions[0])):
        Start_Positions = []
        End_Positions = []
        
        Pos = Seq_Data[Seq_Data["Positions"]==Positions[0][platepositions]]
    
        Start_Pos = [Pos.alignment_genome_start.unique()]
        End_Pos = [Pos.alignment_genome_end.unique()]
        #The following code goes through start positions
        #Teststart gets fixed so it tests how many are at +- 5 bases on each side (Expecting minor errors)
        for copies in range(len(Start_Pos[0])):
            len_total = len(Pos)

            TestStart = Pos[Pos["alignment_genome_start"]>(Start_Pos[0][copies]-5)][Pos["alignment_genome_start"]<(Start_Pos[0][copies]+5)]
            Start_Pos_Predicter = (TestStart["alignment_genome_start"].value_counts())
            Len_Start = len(TestStart["alignment_genome_start"])

            Start_Positions.append([Positions[0][platepositions],Start_Pos_Predicter.index[0],Len_Start/len_total,Pos["Duplicates"].iloc[0]])
        
        
        Start_Positions = np.array([Start_Positions])
        Start_Positions = pd.DataFrame(Start_Positions[0],columns=["PlatePos","Gene_Position","Percentage", "Amount"])
        
        Start_Positions["Percentage"] = pd.to_numeric(Start_Positions["Percentage"])
        #It sorts to put the best choice in the end
        Start_Positions = Start_Positions.sort_values("Percentage")
        ##Here it chooses the best start-position
        Start_Pos_Best = np.array(Start_Positions.iloc[-1,:])

        #Does the same for the end-positions
        for copies in range(len(End_Pos[0])):
            len_total = len(Pos)
            TestEnd = Pos[Pos["alignment_genome_end"]>End_Pos[0][copies]-5][Pos["alignment_genome_end"]<End_Pos[0][copies]+5]
            ###A test to remove the ones that don´t fit with the Fw. primer
            TestEnd = TestEnd[TestEnd["alignment_genome_start"]>int(Start_Pos_Best[1])-5][TestEnd["alignment_genome_start"]<int(Start_Pos_Best[1])+5]
            End_Pos_Predicter = (TestEnd["alignment_genome_end"].value_counts())
            Len_end = len(TestEnd["alignment_genome_end"])
    
            if Len_end>0:
                End_Positions.append([Positions[0][platepositions],End_Pos_Predicter.index[0],Len_end/len_total,Pos["Duplicates"].iloc[0]])
            
            
        
        End_Positions = np.array([End_Positions])
        End_Positions = pd.DataFrame(End_Positions[0],columns=["PlatePos","Gene_Position","Percentage", "Amount"])
        End_Positions["Percentage"] = pd.to_numeric(End_Positions["Percentage"])
        End_Positions = End_Positions.sort_values("Percentage")
        
        End_Pos_Best = np.array(End_Positions.iloc[-1:,:])

        Start_end= np.append(Start_Pos_Best,End_Pos_Best)
        Combined_NGS = np.vstack([Combined_NGS,Start_end])
        
    #This is the full dataframe with start and end-data implemented
    NGS_Full= pd.DataFrame(Combined_NGS.T,["platePos","start_position","start_percentage", "start_amount","platePos","end_position","end_percentage", "end_amount"]).T
    
    #Set data types
    NGS_Full["start_percentage"] =NGS_Full["start_percentage"].astype(float)
    NGS_Full["end_percentage"] =NGS_Full["end_percentage"].astype(float)
    NGS_Full["start_amount"] =NGS_Full["start_amount"].astype(int)
    NGS_Full["end_amount"] =NGS_Full["end_amount"].astype(int)
    NGS_Full["start_position"] =NGS_Full["start_position"].astype(int)
    NGS_Full["end_position"] =NGS_Full["end_position"].astype(int)

    #Filter: Removes the samples which are not found with above parameters. 
    #Might be unnescessary due to the QC in demultiplexer
    NGS_Filtered = NGS_Full[(NGS_Full["start_percentage"]*NGS_Full["start_amount"])>5]
    NGS_Filtered = NGS_Filtered[(NGS_Filtered["end_percentage"]*NGS_Filtered["end_amount"])>5]
    NGS_Filtered = NGS_Filtered[NGS_Filtered["start_percentage"]>0.30]
    NGS_Filtered = NGS_Filtered[NGS_Filtered["end_percentage"]>0.30]
    #Just to make sure that a few point-mutations in the end dont get lost, it is put at orf+5
    #The following removes accidental inserts from the backbone which has been seen rarely. 
    NGS_Filtered = NGS_Filtered[NGS_Filtered["end_position"]<(ORF_Length+5)] 

    
    
    #Luminescence_Implementation############################################################
    #Data from full cell below: 
    
    #Input for this:
    NGS = NGS_Filtered
    Lum =luminescence_path
    Primer = First_primer_used   
    
    FullCell = glob.glob(os.path.join(Lum,"Luciferase_FC_Plate*.xls"))

    Full_Filler = []
    for f in FullCell:
        List_Filter=[] 
        # read the csv file 
        Filter_Data =(pd.read_excel(f).iloc[-96:,2:]).reset_index() #Yields gen5 data

        #A bit of goofy code, should be redone but doesn´t hurt the script speed:
        for i in range(96): 
            List_Filter.append(Filter_Data.iloc[i,1])
        #Redefine positions of data so it fits with "normal" A1, B1...G12, H12 format

        List_Filter2 = []
        for i in range(12):
            for j in range(8):
                List_Filter2.append(List_Filter[(j*12)+(i)])
        #Define plate and loop

        PlateNumber=re.search('_Plate(.*)Time_', f)
        Full_Filler.append([PlateNumber.group(1),List_Filter2])
    #Make data for "Full_Cell"
    Full_Info= []
    Full_columns= []
    
    for i in range(len(Full_Filler)):
        Full_Info.append(Full_Filler[i][1])
        Full_columns.append(Full_Filler[i][0])
    
    Full_Info = np.array(Full_Info).T    
    Full_Info = pd.DataFrame(Full_Info,columns=Full_columns)
    ###End_Full_Cell###

    #Data from cytosol below: 
    Cytosol = glob.glob(os.path.join(Lum,"Luciferase_Cytosol_*.xls"))
    Cytosol_Filler = []
    for f in Cytosol:
        List_Filter=[] 
        # read the csv file 
        Filter_Data =(pd.read_excel(f).iloc[-96:,2:]).reset_index()

        for i in range(96):
            List_Filter.append(Filter_Data.iloc[i,1])
        
        #Redefine positions of data so it fits with "normal" A1, B1...G12, H12 format
        
        List_Filter2 = []
        for i in range(12):
            for j in range(8):
                List_Filter2.append(List_Filter[(j*12)+(i)])
        #Define plate and loop

        
        #Define plate and loop
        PlateNumber=re.search('_Plate(.*)Time_', f)
        Cytosol_Filler.append([PlateNumber.group(1),List_Filter2])
    
    Cytosol_Info= []
    cytosol_columns = []
    #Cytosol_Filler=np.array(Cytosol_Filler)
    for i in range(len(Cytosol_Filler)):
        Cytosol_Info.append(Cytosol_Filler[i][1])
        cytosol_columns.append(Cytosol_Filler[i][0])
    
    Cytosol_Info = np.array(Cytosol_Info).T
    Cytosol_Info = pd.DataFrame(Cytosol_Info,columns=cytosol_columns)
    ###End_Cytosol###
    #Combine above to a readable dataframe
    
    Insertion_array = np.array([0,0,0])
    for i in range(len(Cytosol)):
        Cytosol_plate = Cytosol_Info[f'{i+1}']
        Full_plate = Full_Info[f'{i+1}']
        for j in range(96):
            Position =  f"{i+Primer}"+","+f"{j+1}"
            Cytosol_Lum = Cytosol_plate[j]
            Full_Lum = Full_plate[j]
            Insertion_point = np.array([Position,Cytosol_Lum,Full_Lum])
            Insertion_array = np.vstack([Insertion_array,Insertion_point])

    Insertion_array = pd.DataFrame(Insertion_array,columns=["platePos","Cytosol","Full_Cell"])
    Insertion_DF=Insertion_array[Insertion_array["platePos"]!=0]

    NGS_Filtered = NGS.iloc[:,1:]
    NGS_Filtered["platePos"] = NGS_Filtered["platePos"].astype(str)
    Insertion_DF["platePos"] = Insertion_DF["platePos"].astype(str)

    NGS_Filtered = NGS_Filtered.merge(Insertion_DF, left_on='platePos', right_on ='platePos' )

    ###############This part extracts OD-data and implements it in the data#############
    pattern = "OD_Plate*.xlsx"
    # Use glob to get all the files matching the pattern
    files = glob.glob(f"{luminescence_path}\\{pattern}")


    ####Data from growth cultures###########
    # Create a dictionary to store all the dataframes with their coresponding file names as keys
    dfs = {}

    # Read each excel file into a pandas DataFrame and store it in the dictionary
    for file in files:
        # Extract the file number from the file name to use as a key
        file_num = file.split('OD_Plate')[1].split('.xlsx')[0]
        
        # Read the Excel file into a DataFrame
        dfs[file_num] = pd.read_excel(file)


    letters = ["A", "B", "C", "D", "E", "F", "G", "H"]

    OD_data = np.array(["","","","",""])

    for plates in range(len(dfs)):
        posnumber = 1
        for j in range(12):
            for i in range(8):
                #Get the time and 
                A= dfs[f"{plates+1}"]["Time"][:5]
                B= dfs[f"{plates+1}"][f"{letters[i]}{j+1}"][:5]
                plt.plot(A,B)
                trap1=np.trapz(B,A)
                #print(trap)
                
                A= dfs[f"{plates+1}"]["Time"][4:]
                B= dfs[f"{plates+1}"][f"{letters[i]}{j+1}"][4:]
                plt.plot(A,B)
                trap2=np.trapz(B-dfs[f"{plates+1}"][f"{letters[i]}{j+1}"][4],A)
                
                OD_data = np.vstack((OD_data, np.array([f"{letters[i]}{j+1}",plates+Primer,f"{plates+Primer,posnumber}",trap1,trap2])))
                
                posnumber+=1
            
    #############Here is pre-OD-Data

    #Data from pre-culture:
    Pre_exp = glob.glob(f"{luminescence_path}\\Pre_Experiment_Data_Loop__Plate_*.xls")
    pre_temp = []
    predata_total = np.array(["","",""])
    for f in Pre_exp:
        List_Filter=[] 
        List_pos = []
        # read the csv file
        predata_plate =pd.read_excel(f).iloc[-96:,3]
        if len(predata_plate)==96:
            for i in range(96):
                List_Filter.append(predata_plate.iloc[i])
        
        predata_pos =pd.read_excel(f).iloc[-96:,2]
        if len(predata_pos)==96:
            for i in range(96):
                List_pos.append(predata_pos.iloc[i])

            #Define plate and loop
            Loopnumber=re.search('_Loop_(.*)_', f)
            ###Adding "Primer-1" to fix if the second primerset is used
            PlateNumber= re.search('_Plate_(.*)Time_Date_', f)

            Time_Date = re.search('Time_Date_(.*).xls', f)
            #Extract data from above and combine it (Primer is added to fit with "Real" position)
            pre_temp=([int(PlateNumber.group(1))+Primer-1,List_pos,List_Filter])

            predata=pd.DataFrame(np.array([np.repeat(pre_temp[0],96),pre_temp[1],pre_temp[2]])).T
        predata_total = np.vstack((predata_total,predata))

    od_df = pd.DataFrame(OD_data,columns=["pos","plate","platePos","preinduction_auc","postinduction_auc"])
    pre_df = pd.DataFrame(predata_total,columns=["plate","pos","pre_od"])
    merged_od = pd.merge(od_df, pre_df, on=['plate', 'pos'])
    ###Prep the platePos  for merge
    merged_od['platePos'] = merged_od['platePos'].str.replace('(', '').str.replace(')', '').str.replace(' ','')

    NGS_Filtered = pd.merge(merged_od, NGS_Filtered, on=['platePos'])    
    NGS_Filtered.to_csv(f'{Output_Path}/NGS_Filtered.csv', index = False)

    return NGS_Filtered
  
##Problem: Loads twice on calm?### test random stuff in "return"