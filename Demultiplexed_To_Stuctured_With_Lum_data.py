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
        #print(List_Filter2)

        PlateNumber=re.search('_Plate(.*)Time_', f)
        Full_Filler.append([PlateNumber.group(1)[0],List_Filter2])
    #Make data for "Full_Cell"
    Full_Info= []
    
    Full_Filler=np.array(Full_Filler)
    for i in range(16):
        Full_Info.append(Full_Filler[i][1])
    
    Full_Info = np.array(Full_Info).T    
    Full_Info = pd.DataFrame(Full_Info,columns=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"])
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
        Cytosol_Filler.append([PlateNumber.group(1)[0],List_Filter2])
    
    Cytosol_Info= []
    
    Cytosol_Filler=np.array(Cytosol_Filler)
    for i in range(16):
        Cytosol_Info.append(Cytosol_Filler[i][1])
    
    Cytosol_Info = np.array(Cytosol_Info).T
    Cytosol_Info = pd.DataFrame(Cytosol_Info,columns=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"])
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

    ######################################
    
    #Data from pre-culture:

    Pre_exp = glob.glob(os.path.join(luminescence_path,"Pre_Experiment_Data_Loop__Plate_*.xls"))

    PreData = []
    for f in Pre_exp:
        List_Filter=[] 
        # read the csv file
        #Error: Data kommer ud 2 før i L2P2
        Filter_Data =pd.read_excel(f).iloc[-96:,3]
        if len(Filter_Data)==96:
            for i in range(96):
                List_Filter.append(Filter_Data.iloc[i])
    
    
            #Define plate and loop
            Loopnumber=re.search('_Loop_(.*)_', f)
            PlateNumber=re.search('_Plate_(.*)Time_Date_', f)
    
            Time_Date = re.search('Time_Date_(.*).xls', f)
    
            PreData.append([Time_Date.group(1)[:],PlateNumber.group(1)[0],Loopnumber.group(1)[0],List_Filter])
            
    
    
    Pre_ODs = pd.DataFrame([])
    positions = []
    for i in range(8):
        for j in range(12):
            positions.append((1+i)+(j*8))
    
    positions=np.array(positions)
    
    for i in range(16):
        df = pd.DataFrame(PreData[i][3],columns=["preODs"])
        df["plate"]=i+1
        df["position"]= positions
        Pre_ODs = Pre_ODs.append(df)
    
    Pre_ODs["platePos"] = Pre_ODs['plate'].astype(str)+","+Pre_ODs['position'].astype(str)

    ####Now for the plotting!
    
    NGS=NGS_Filtered.sort_values("Cytosol")
    
    #print(NGS_Filtered.head())
    #Plotter: 
    
    #Go through the data and extract all the start/end/size data: 
    Data_Distribution = []
    for i in range((NGS["start_percentage"].size)):
        Data_Distribution.append([NGS_Filtered["start_position"].iloc[i],NGS_Filtered["end_position"].iloc[i],NGS_Filtered["end_position"].iloc[i]-NGS_Filtered["start_position"].iloc[i]])
    #Create the data for plotting:
    Plotting_info=[]
    #Make Data_Distribution into dataframe
    Data_Distribution= np.asarray(Data_Distribution)
    
    Data_Distribution=pd.DataFrame(Data_Distribution,columns=["Start","Stop","Length"])
    
    #Go through the size of Data_Distribution, At every point gain size-data for each base
    for i in range(len(Data_Distribution["Start"])):
        for j in range(Data_Distribution["Start"][i],Data_Distribution["Stop"][i]):
            Full_cell=np.array([NGS_Filtered["Full_Cell"][i]]).astype(float)
            Cytosol=np.array([NGS_Filtered["Cytosol"][i]]).astype(float)
            Plotting_info.append([j,Data_Distribution["Length"][i],Full_cell,Cytosol])
            #break
        #break
    #Make plotting info into dataframe
    Plotting_info2=np.asarray(Plotting_info)
    Plotting_info2=pd.DataFrame(Plotting_info2,columns=["Base","Length","Full_Cell","Cytosol"]).sort_values("Cytosol")
    
    #Extra high outlier removal below: 
    #Plotting_info2=Plotting_info2[Plotting_info2["Full_Cell"]<1000]
    #Plotting_info2=Plotting_info2[Plotting_info2["Cytosol"]<1000]
    
    
    #Colour the plot by luminescence
    
    #Plot the thing! 
    color_norm = mcolors.LogNorm(vmin=Plotting_info2["Full_Cell"].min(), vmax=Plotting_info2["Full_Cell"].max())
    plt.scatter(Plotting_info2["Base"],Plotting_info2["Length"],c=Plotting_info2["Full_Cell"], cmap=plt.cm.get_cmap('coolwarm', 10),s=0.1,norm=color_norm)
    plt.title("Full cell_Lum")
    plt.xlabel("Base_Number")
    plt.ylabel("Size")
    plt.xlim([0,ORF_Length])
    plt.ylim([0,ORF_Length])
    plt.colorbar()

    plt.savefig(f"{Output_Path}/FC_Lum.png",dpi=600)
    plt.show()
    
    color_norm = mcolors.LogNorm(vmin=Plotting_info2["Cytosol"].min(), vmax=Plotting_info2["Cytosol"].max())
    plt.scatter(Plotting_info2["Base"],Plotting_info2["Length"],c=Plotting_info2["Cytosol"], cmap=plt.cm.get_cmap('coolwarm', 10),s=0.1,norm=color_norm)
    plt.title("Cytosol_Lum")
    plt.xlabel("Base_Number")
    plt.ylabel("Size")
    plt.xlim([0,ORF_Length])
    plt.ylim([0,ORF_Length])
    
    plt.colorbar()

    plt.savefig(f"{Output_Path}/Cytosol_lum.png",dpi=600)
    plt.show()

    #Data from above but times the length to get an approximate amount of protein: 
    NGS_Filtered["FC_Size"] = NGS_Filtered["Full_Cell"].astype(int)*(NGS_Filtered["end_position"]-NGS_Filtered["start_position"])/10000
    
    NGS_Filtered["C_Size"] = NGS_Filtered["Cytosol"].astype(int)*(NGS_Filtered["end_position"]-NGS_Filtered["start_position"])/10000
    
    Data_Distribution = []
    for i in range((NGS["start_percentage"].size)):
        Data_Distribution.append([NGS_Filtered["start_position"].iloc[i],NGS_Filtered["end_position"].iloc[i],NGS_Filtered["end_position"].iloc[i]-NGS_Filtered["start_position"].iloc[i]])
    #Create the data for plotting:
    Plotting_info=[]
    #Make Data_Distribution into dataframe
    Data_Distribution= np.asarray(Data_Distribution)
    
    Data_Distribution=pd.DataFrame(Data_Distribution,columns=["Start","Stop","Length"])
    
    #Go through the size of Data_Distribution, At every point gain size-data for each base
    for i in range(len(Data_Distribution["Start"])):
        for j in range(Data_Distribution["Start"][i],Data_Distribution["Stop"][i]):
            Full_cell=((np.array([NGS_Filtered["FC_Size"][i]]).astype(float)))
            Cytosol=((np.array([NGS_Filtered["C_Size"][i]]).astype(float)))
            Plotting_info.append([j,Data_Distribution["Length"][i],Full_cell,Cytosol])
        #if i>200: #Break to not break my computer :D 
        #    break
    #Make plotting info into dataframe
    Plotting_info2=np.asarray(Plotting_info)
    Plotting_info2=pd.DataFrame(Plotting_info2,columns=["Base","Length","Full_Cell","Cytosol"])
    
    #Go through the size of Data_Distribution, At every point gain size-data for each base
    for i in range(len(Data_Distribution["Start"])):
        for j in range(Data_Distribution["Start"][i],Data_Distribution["Stop"][i]):
            Full_cell=((np.array([NGS_Filtered["FC_Size"][i]]).astype(float)))
            Cytosol=((np.array([NGS_Filtered["C_Size"][i]]).astype(float)))
            Plotting_info.append([j,Data_Distribution["Length"][i],Full_cell,Cytosol])

    #Make plotting info into dataframe
    Plotting_info2=np.asarray(Plotting_info)
    Plotting_info2=pd.DataFrame(Plotting_info2,columns=["Base","Length","Full_Cell","Cytosol"]).sort_values("Cytosol")
    

    
    #Colour the plot by luminescence
    
    #Plot the thing! 
    color_norm = mcolors.LogNorm(vmin=(Plotting_info2["Full_Cell"]*Plotting_info2["Length"]/10000).min(), vmax=(Plotting_info2["Full_Cell"]*Plotting_info2["Length"]/10000).max())
    plt.scatter(Plotting_info2["Base"],Plotting_info2["Length"],c=Plotting_info2["Full_Cell"]*Plotting_info2["Length"]/10000, cmap=plt.cm.get_cmap('coolwarm', 10),s=0.01)
    plt.title("Full cell_LumXsize")
    plt.xlabel("Base_Number")
    plt.ylabel("Size")
    plt.xlim([0,ORF_Length])
    plt.ylim([0,ORF_Length])
    plt.colorbar()
    #plt.clim(0,10)
    plt.savefig(f"{Output_Path}/FC_LumXsize.png",dpi=600)
    plt.show()
    
    
    color_norm = mcolors.LogNorm(vmin=(Plotting_info2["Cytosol"]*Plotting_info2["Length"]/10000).min(), vmax=(Plotting_info2["Cytosol"]*Plotting_info2["Length"]/10000).max())
    plt.scatter(Plotting_info2["Base"],Plotting_info2["Length"],c=Plotting_info2["Cytosol"]*(Plotting_info2["Length"])/10000, cmap=plt.cm.get_cmap('coolwarm', 10),s=0.01)
    plt.title("Cytosol_LumXsize")
    plt.xlabel("Base_Number")
    plt.ylabel("Size")
    plt.xlim([0,ORF_Length])
    plt.ylim([0,ORF_Length])
    plt.colorbar()
    #plt.clim(0,5)
    plt.savefig(f"{Output_Path}/C_LumXsize.png",dpi=600)
    plt.show()
    ###END of input
    
    
    NGS_Filtered.to_csv(f'{Output_Path}/NGS_Filtered.csv', index = False)
    NGS_Full.to_csv(f'{Output_Path}/NGS_Full.csv', index = False)
    
    NGS_Filtered["Cytosol"]=NGS_Filtered["Cytosol"].astype(int)
    NGS_Filtered["Full_Cell"]=NGS_Filtered["Full_Cell"].astype(int)
    NGS_Filtered=NGS_Filtered.sort_values("Cytosol")
    
    ####Data to make freq-plots:
    combined_fc= np.zeros(ORF_Length)
    combined_c= np.zeros(ORF_Length)
    amount_seq = np.zeros(ORF_Length)
    for i in range(len(NGS_Filtered["start_position"])):
        combined_fc[NGS_Filtered["start_position"][i]:NGS_Filtered["end_position"][i]] += NGS_Filtered["Full_Cell"][i]
        combined_c[NGS_Filtered["start_position"][i]:NGS_Filtered["end_position"][i]] += NGS_Filtered["Cytosol"][i]
        amount_seq[NGS_Filtered["start_position"][i]:NGS_Filtered["end_position"][i]]+=1

        
    xAxis = np.arange(0,ORF_Length)

    plt.plot(xAxis,combined_fc/amount_seq, linewidth=0.2)
    plt.plot(xAxis,combined_c/amount_seq, linewidth=0.2)
    plt.scatter(xAxis,combined_fc/amount_seq,c=amount_seq, cmap='coolwarm',s=1)
    plt.scatter(xAxis,combined_c/amount_seq,c=amount_seq, cmap='coolwarm',s=1)
    #print(len(amount_seq),len(combined_c))

    
    plt.title("FC_C_prBase")
    plt.xlabel("Base_Number")
    plt.ylabel("Luminescence_x_base")
    plt.colorbar()
    #plt.clim(0,100)
    plt.savefig(f"{Output_Path}/Base_Composition.png",dpi=600)
    
    plt.show()

    
    combined_pos= np.zeros(ORF_Length)

    for i in range(len(NGS_Filtered["start_position"])):
        combined_pos[NGS_Filtered["start_position"][i]:NGS_Filtered["end_position"][i]] += (NGS_Filtered["C_Size"][i])/(NGS_Filtered["FC_Size"][i])

      
        
    

    plt.scatter(xAxis,combined_pos/amount_seq,c=amount_seq, cmap='coolwarm',s=1)
    plt.scatter(xAxis,combined_fc/amount_seq,c=amount_seq, cmap='coolwarm',s=1)
    plt.scatter(xAxis,combined_c/amount_seq,c=amount_seq, cmap='coolwarm',s=1)
    plt.title("C/FC")
    plt.xlabel("Base_Number")
    plt.ylabel("Approx_Fraction")
    plt.colorbar()
    #plt.clim(0,100)
    plt.savefig(f"{Output_Path}/solubility.png",dpi=600)
    
    plt.show()


        # Count barcodes
    filtered = NGS_Filtered
    start_amounts = np.array(filtered["start_amount"])
    normed_amounts = (start_amounts - start_amounts.min()) / (start_amounts.max() - start_amounts.min()) # Normalize values between 0 and 1
    colormap = cm.get_cmap('coolwarm') # Choose a colormap; you can use others like 'plasma', 'inferno', etc.
    
    for i in range(len(filtered["start_amount"])):
        color = colormap(normed_amounts[i])
        plt.scatter(int(filtered["platePos"][i].split(',')[1]),int(filtered["platePos"][i].split(',')[0]), c=[color], s=25)
    
    
    plt.xlabel("Fw_primers")
    plt.ylabel("Rev_primers")
    
    colorbar = plt.colorbar(cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(start_amounts.min(), start_amounts.max())))
    colorbar.set_label('Reads found')
    
    plt.savefig(f"{Output_Path}/barcoding_Description", dpi=600)
    plt.show()

    
    return NGS_Filtered

