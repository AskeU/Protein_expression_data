# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:47:01 2022

@author: askung
"""
##PACKAGES##
import pandas as pd
import numpy as np
import os
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns

##Find path##
path = "C:/Users/askung/OneDrive - Danmarks Tekniske Universitet/WP1/Automation_Work/2023_11_16_p85a_1000rpm_37"

plt.rcParams['font.size'] = 25

#Find relevant files for OD and GFP
Pre_Induction = glob.glob(os.path.join(path,"Pre_Induction_Data_Loop_*.xls"))
Induction_Data = glob.glob(os.path.join(path,"Induction*.xls"))

#Set WD
os.chdir(path)


#Data from pre-culture:
Platenumbers=[]
Loopnumbers=[]
PreData = []
for f in Pre_Induction:
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

        PreData.append([Time_Date.group(1)[:],PlateNumber.group(1)[0:2],Loopnumber.group(1)[0],List_Filter])



#Data from induction
Platenumbers=[]
Loopnumbers=[]
Induction = []
for f in Induction_Data:
    List_Filter=[] 
    # read the csv file
    Filter_Data =pd.read_excel(f).iloc[-96:,3:]
    if len(Filter_Data)==96:
        for i in range(96):
            List_Filter.append(Filter_Data.iloc[i,0])
    
    
        #Define plate and loop
        Loopnumber=re.search('_Loop_(.*)_', f)
        PlateNumber=re.search('_Plate_(.*)Time_Date_', f)
        Time_Date = re.search('Time_Date_(.*).xls', f)
        Induction.append([Time_Date.group(1)[:],PlateNumber.group(1)[0:2],Loopnumber.group(1)[0],List_Filter])
    
PreData = pd.DataFrame(PreData,columns=["Time","Plate","Loop","Data"])
InductionData = pd.DataFrame(Induction,columns=["Time","Plate","Loop","Data"])

Full_Data = pd.concat([PreData,InductionData]).reset_index()
Full_Data =Full_Data.drop(['index'], axis=1)

#Get the start-timer

Startmin = int(PreData["Time"][0][3:5])
Starthour = int(PreData["Time"][0][0:2])
Startday = int(PreData["Time"][0][-2:])
Startmonth = int(PreData["Time"][0][-5:-3])
Startyear = int(PreData["Time"][0][-8:-6])

for i in range(len(Full_Data["Time"])):
    test_min = (int(Full_Data["Time"][i][3:5])-Startmin)/60
    test_hour = int(Full_Data["Time"][i][0:2])-Starthour
    test_day = (int(Full_Data["Time"][i][-2:])-Startday)*24
    test_month = (int(Full_Data["Time"][i][-5:-3])-Startmonth)
    test_year = int(Full_Data["Time"][i][-8:-6])-Startyear
    
    #Doesn´t work if experiment runs multiple days over end-month
    if test_month<0:
        test_day = test_day+Startday
    

    
    Full_Data["Time"][i]=test_min+test_hour+test_day
    

Columns=[]
Letters = ["A","B","C","D","E","F","G","H"]
for i in Letters:
    for j in range(12):
        Columns.append(f'{i}{j+1}')
Columns.append("Time")

#Data from full cell below: 

FullCell = glob.glob(os.path.join(path,"Luciferase_FC_Plate*.xls"))
Full_Filler = []
for f in FullCell:
    List_Filter=[] 
    # read the csv file 
    Filter_Data =(pd.read_excel(f).iloc[-96:,2:]).reset_index()

    for i in range(96):
        List_Filter.append(Filter_Data.iloc[i,1])

    #Define plate and loop

    PlateNumber=re.search('_Plate(.*)Time_', f)
    Full_Filler.append([PlateNumber.group(1)[0],List_Filter])
#Make data for "Full_Cell"
Full_Info= []

#Full_Filler=np.array(Full_Filler)

for i in range(len(Full_Filler)):
    Full_Info.append(Full_Filler[i][1])

Full_Info = np.array(Full_Info)    
###End_Full_Cell###

#Data from cytosol below: 

Cytosol = glob.glob(os.path.join(path,"Luciferase_Cytosol_*.xls"))
Cytosol_Filler = []
for f in Cytosol:
    List_Filter=[] 
    # read the csv file 
    Filter_Data =(pd.read_excel(f).iloc[-96:,2:]).reset_index()

    for i in range(96):
        List_Filter.append(Filter_Data.iloc[i,1])

    #Define plate and loop

    PlateNumber=re.search('_Plate(.*)Time_', f)
    Cytosol_Filler.append([PlateNumber.group(1)[0],List_Filter])

Cytosol_Info= []

#Cytosol_Filler=np.array(Cytosol_Filler)
for i in range(len(Cytosol_Filler)):
    Cytosol_Info.append(Cytosol_Filler[i][1])

Cytosol_Info = np.array(Cytosol_Info)    
###End_Cytosol###

##Start TrxA_control##
TrxA_Control = glob.glob(os.path.join(path,"TrxA_Control_Time*.xls"))
Filter_Data =(pd.read_excel(TrxA_Control[0]).iloc[31:,2:]).reset_index()
Control_Data =np.array(Filter_Data[9:17])

#Make plots of the data, start with growth curves
sns.set(rc = {'figure.figsize':(20,10)})
for plates in range(int(Full_Data["Plate"].astype(int).max())):
    sns.set(font_scale = 2)
    Single_Plate = (pd.DataFrame(Full_Data[Full_Data["Plate"]==str(plates+1)])).reset_index()
    XAxis = pd.DataFrame(Single_Plate["Time"])
    
    
    YAxisFull = []
    for j in range(96):
        YAxis_Pos_Collecter = []
        for i in range(len(Single_Plate["Time"])):
            YAxis_Pos_Collecter.append(Single_Plate["Data"][i][j])
        YAxisFull.append(YAxis_Pos_Collecter)



    Export_File = (pd.DataFrame(YAxisFull).T)
    Export_File["Time"] = XAxis

    Export_File=Export_File.sort_values("Time")
    Export_File.columns=Columns
    Export_File.to_excel(f'OD_Plate{plates+1}.xlsx')
    

    for positions in range(96):
        plt.plot(Export_File["Time"],Export_File[f'{Columns[positions]}']) # plotting t, a separately 
    plt.ylim(ymax = 1.8, ymin = 0)    
    plt.title(f'Plate{plates+1}')
    plt.show()
    
###############Get average of last 3 OD reads##############
    p1,p2 =-3,-1
    Heatmaps = []
    df= []    
    Index= ["A","B","C","D","E","F","G","H"]
    Cols = []
    
    Export_File=Export_File.reset_index()
    
    for columns in range(12):
        Column = (columns)*8
        df.append([Export_File[f'A{columns+1}'][p1:p2].mean(),Export_File[f'B{columns+1}'][p1:p2].mean(),Export_File[f'C{columns+1}'][p1:p2].mean(),Export_File[f'D{columns+1}'][p1:p2].mean(),Export_File[f'E{columns+1}'][p1:p2].mean(),Export_File[f'F{columns+1}'][p1:p2].mean(),Export_File[f'G{columns+1}'][p1:p2].mean(),Export_File[f'H{columns+1}'][p1:p2].mean()])
    
    df=pd.DataFrame(df).T
    


    sns.heatmap(df, vmin=0, vmax=2, annot=True,  cmap="RdYlGn" ,yticklabels= Index,xticklabels= [1,2,3,4,5,6,7,8,9,10,11,12]).set(title=f'End_OD: Plate{plates+1}')
    plt.show()



########Also see OD about middle of experiment###############

    p1,p2 =10,12
    Heatmaps = []
    df= []    
    Index= ["A","B","C","D","E","F","G","H"]
    Cols = []
    
    Export_File=Export_File.reset_index()
    
    for columns in range(12):
        Column = (columns)*8
        df.append([Export_File[f'A{columns+1}'][p1:p2].mean(),Export_File[f'B{columns+1}'][p1:p2].mean(),Export_File[f'C{columns+1}'][p1:p2].mean(),Export_File[f'D{columns+1}'][p1:p2].mean(),Export_File[f'E{columns+1}'][p1:p2].mean(),Export_File[f'F{columns+1}'][p1:p2].mean(),Export_File[f'G{columns+1}'][p1:p2].mean(),Export_File[f'H{columns+1}'][p1:p2].mean()])
    
    df=pd.DataFrame(df).T

    

    sns.heatmap(df, vmin=0, vmax=2, annot=True, cmap="RdYlGn", yticklabels=Index, xticklabels=[1,2,3,4,5,6,7,8,9,10,11,12]).set_title(f'Mid_OD: Plate{plates+1}', fontsize=25)


    plt.show()


    #Make heatmap of full_Cell
    Number_Pos =0
    Heat_FC = np.zeros(shape=(8,12))
    for rows in range(8):
        for columns in range(12):
            Heat_FC[rows,columns]=  Full_Info[plates][Number_Pos]
            Number_Pos=Number_Pos+1
            df=pd.DataFrame(Heat_FC)

    sns.set(font_scale = 1.3)
    
    sns.heatmap(df, vmax=np.max(Control_Data), annot=True, cmap="RdYlGn", yticklabels=Index, xticklabels=[1,2,3,4,5,6,7,8,9,10,11,12]).set_title(f'Full_Lum: Plate{plates+1}', fontsize=25)

    plt.show()
#Make heatmap of cytosol
    Number_Pos =0
    Heat_C = np.zeros(shape=(8,12))
    for rows in range(8):
        for columns in range(12):
            Heat_C[rows,columns]=  Cytosol_Info[plates][Number_Pos]
            Number_Pos=Number_Pos+1
            df=pd.DataFrame(Heat_C)

    
    
    sns.heatmap(df, vmax=np.max(Control_Data), annot=True, cmap="RdYlGn", yticklabels=Index, xticklabels=[1,2,3,4,5,6,7,8,9,10,11,12]).set_title(f'Cytosol_Lum: Plate{plates+1}', fontsize=25)

    plt.show()
#Heatmap of solubility level
    Solubility = Heat_C/Heat_FC
    Solubility[Heat_FC<100]= -1
    

    sns.set(font_scale = 2)
    
    sns.heatmap(Solubility, vmin=0, vmax=1, annot=True, cmap="RdYlGn", yticklabels=Index, xticklabels=[1,2,3,4,5,6,7,8,9,10,11,12]).set_title(f'Solubility: Plate{plates+1}', fontsize=25)

    plt.show()

Testing =np.array([Full_Info.ravel(),Cytosol_Info.ravel()]).T
TestingPD = (pd.DataFrame(Testing, columns = ['Full_Cell','Cytosolic'])).sort_values(by = "Full_Cell")

#Testing = np.sort(Full_Info.ravel())
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (10,10))
Xvalue = np.arange(1,TestingPD.shape[0]+1)
sns.scatterplot(Xvalue,TestingPD['Full_Cell'], ax=ax, color='C0', label = 'Full_Cell')
sns.scatterplot(Xvalue,TestingPD['Cytosolic'], ax=ax, color='C1', label = 'Cytosolic')
plt.legend()

plt.show()


