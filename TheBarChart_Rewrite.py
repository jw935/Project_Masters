import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from decimal import Decimal as D

Gconsts = {'DEBUG' : False,
           'PixConv' : 1/1280, # To convert pixels to metres
           'FontSize' : 12, # Size of font on the graph generated
           'VertNum' : 0, # Number assigned to each disc starts at 0 at S1-L5 disc then counts up to 8
           'Save' : False, # Save the graph?
           'Labels' : ['L5','L4','L3','L2','L1','T12','T11','T10'], # Vertebral names
           'NumofVert' : 9,
           'time' : 10,
           'E1' : 90.42*10**9,
           'c' : -1553.26*10**9,
           'Area' : [0.0008356799406036492,0.0008783392577201085,0.0008686146838902883,0.000876787817241108,
                   0.000794520980955232,0.0007278585345006305,0.0007038109103025612,0.0006719719699558038]
           # Areas of the discs average in m^2. All these consts coppied from outputs of previous scripts.
           }   

if Gconsts['DEBUG'] == True:
    print("Current working directory: ", format(os.getcwd()))

os.chdir('I:\MPhys Project\Data')
FileNs = np.genfromtxt('testfiles.txt', dtype=str)

if Gconsts['DEBUG'] == True:
    print(FileNs)
     
os.chdir('I:\MPhys Project\Data\Centres_of_Rotation')

if Gconsts['DEBUG'] == True:
    print("Current working directory: ", format(os.getcwd()))

COR8 = np.genfromtxt('8kgCoR.txt')
COR16 = np.genfromtxt('16kgCoR.txt')

os.chdir('I:\MPhys Project\Data\LoadedSpineData')

if Gconsts['DEBUG'] == True:
    print("Current working directory: ", format(os.getcwd()))
    print('-------------------------------------')

def ReadIn(Data, Gconsts):
    Horr=[]
    Vert=[]
    if Gconsts['DEBUG'] == True:
        print(Data)
        print("----------")
        
    # Convertes the values from pixels to metres and inverts the y axis. Vert is the y componant of the points Horr is the x componant. 
    # These are plotted with Vert on the x axis and Horr on the y axis due to the monge parameterisation.
    for i in range(len(Data)):
        Vert.append(float(Data[i][0])*Gconsts['PixConv'])
        Horr.append(float(Data[i][1])*Gconsts['PixConv'])
        
    for i in range (len(Horr)):
        Horr[i]=-1*Horr[i]
        
    vertCentreBot=str(Vert[14]) # 14th point in each set is the centre of the bottom endplate of the L5 disc.  
    horrCentreBot=str(Horr[14]) # This is moved to (0,0) for every dataset.
    
    if Gconsts['DEBUG'] == True:
        print(horrCentreBot,vertCentreBot)
        print(Horr,Vert)
        print("----------")
        
    for i in range (len(Vert)):
        # Makes the centre of the sacrum endplate (Point 14) the zero point and deals with the floating point error.
        Vert[i]=float(D(str(Vert[i]))-D(vertCentreBot))
        Horr[i]=float(D(str(Horr[i]))-D(horrCentreBot))
        
    if Gconsts['DEBUG'] == True:
        print(Horr)
        print(Vert)
        print("----------")
    return Horr,Vert

def Centering(File, Horr, Vert):
    # Gives a list of x's and y's for the center points of each endplate 
    cHorr=[]
    cVert=[]
    # c prefix for centered
    for i in range(int(len(Horr)/14)):
        cHorr.append(Horr[i*14])
        cVert.append(Vert[i*14])
        
    if Gconsts['DEBUG'] == True:
        print("Centre Coords")
        print(cHorr,cVert)   
    return cHorr,cVert 

def Endplate_Centering(File, Horr, Vert):
    # Brings the centered coords together and puts them in a dictionary for easier access.
    cHorr,cVert = Centering (File, Horr, Vert)
    keys = Gconsts['Labels']
    Centers = []
    for i in range (Gconsts['NumofVert']):
        Centers.append([(cHorr[2*i],cVert[2*i]),(cHorr[2*i+1],cVert[2*i+1])])
        
    if Gconsts['DEBUG'] == True:
        print(Centers)  
        
    Endplate_Centers = dict(zip(keys,Centers))
    if Gconsts['DEBUG'] == True:
        print(Endplate_Centers)
        print('-------------------------------------')
    return Endplate_Centers,keys

def Endplate_Points(Horr,Vert,Vertebra): 
    # Returns two dictionaries of the points accross the top and bootom endplates of a given IVD.
    PointNum = np.arange(0,9,1) # The number assigned to each point along the upper and lower endplates
    TopCoords = []              # so to match corrosponding points and identify in list.
    BottomCoords = []
    for i in range (9):
        BottomCoords.append((Horr[10+i],Vert[10+i]))
    n = 0
    for i in range (9):
        if n == 5: # Required as the 0th endplate point is in the middle of the endplate in original files.
            n = n - 28  
        TopCoords.append((Horr[32-n],Vert[32-n]))
        n=n+1
        
    Top = dict(zip(PointNum,TopCoords))
    Bottom = dict(zip(PointNum, BottomCoords))
    if Gconsts['DEBUG'] == True:
        print(Vertebra)
        print(Top)
        print()
        print(Bottom)
        print('-------------------------------------')
        
    return Top, Bottom # Coords along the top and bottom endplates of an IVD respectively. (Dict)

def Pythag(x,y):
    #Just a Pythagoras func.
    r=((x**2 + y**2)**0.5)
    return r

def Length_Adverager(top,bottom):
    #Outputs a list of all the widths accross the disc between the corosponding points.
    l_List = []
    for i in range(len(top)):
        deltax = float(D(str(top[i][0])) - D(str(bottom[i][0])))
        deltay = float(D(str(top[i][1])) - D(str(bottom[i][1])))
        Dist = Pythag(deltax, deltay)
        l_List.append(Dist)
        lengthAvg = np.average(l_List)
    return lengthAvg

def Length_Finder(Vertebra,File):
    #Finds L and Lo which are the loaded and unloaded lengths of the IVD respectively 
    #by taking an adverage of the output of length_Adverager.
    Data0 = np.genfromtxt(File, dtype=str, skip_header=3, skip_footer=1)
    Horr,Vert = ReadIn(Data0, Gconsts)
    Top, Bottom = Endplate_Points(Horr,Vert,Vertebra)
    l = Length_Adverager(Top, Bottom)
    lmm = (l)*10**3
    return lmm

'''ACTUAL VALUE'''

L8_list = []
L16_list = []

if Gconsts['DEBUG'] == True:
    print('--',Gconsts['Labels'][Gconsts['VertNum']],'Actual','--')

for SubjectNum in range(5): # Cycling the 5 participants in the test set.
    Participant = 3*SubjectNum
    Dataset = FileNs[Participant]
    Dataset8 = FileNs[Participant + 1]
    Dataset16 = FileNs[Participant + 2]    
    Ll8 = Length_Finder(Gconsts['Labels'], Dataset8)
    Ll16 = Length_Finder(Gconsts['Labels'], Dataset16)
    
    L8_list.append(Ll8)
    L16_list.append(Ll16)
    
    if Gconsts['DEBUG'] == True:
        print('--',FileNs[Participant],'--') 
        print('Ll 8kg (mm) = ',Ll8)
        print('Ll 16kg (mm) = ', Ll16)

eao, eat, eath, eaf, eafi = L8_list  # (e.g. eafi = 'E'ight kg 'A'cctual test set 'FI've)
sao, sat, sath, saf, safi = L16_list


'''PREDICTED VALUE'''

Lo_Bar = []

def predicted_strain(stress,Gconsts): # finding predicted strain using equation derrived elsewhere (SLSM)
    strain = (stress/Gconsts['E1'])* (1-np.exp(-1*Gconsts['E1']*Gconsts['time']/Gconsts['c']))
    return strain

for n in range(2):
    # Loop over the two loads. Load number = n + 1
    stress = (n + 1)*8*9.81/Gconsts['Area'][Gconsts['VertNum']]
    if Gconsts['DEBUG'] == True:
        print('--',Gconsts['Labels'][Gconsts['VertNum']],(n + 1)*8,'kg Predicted','--')
        
    Lload_list = []
    
    for k in range(5):
        SubjectNum = k
        Twit = 3*SubjectNum
        Dataset = FileNs[Twit]
        Dataset8 = FileNs[Twit + 1]
        Dataset16 = FileNs[Twit + 2]
        # obtaining initial length
        Lo = Length_Finder(Gconsts['Labels'][Gconsts['VertNum']], Dataset)
        Pstrain = predicted_strain(stress,Gconsts)
        dL = Pstrain/Lo # dl = change in L
        LFinal = dL + Lo
        if Gconsts['DEBUG'] == True:
            print('dL', dL)
        Lload_list.append(LFinal)
        if n == 1:
            Lo_Bar.append(Lo)
            Lo_Bar.append(Lo) # Done twice to obtain a green line that spans 2 groups of bars denoting the original length of the IVD.

    if (n + 1) == 1:
        epo, ept, epth, epf, epfi = Lload_list # Two to catch the two loading scenarios (8kg or 16kg).
    if (n + 1) == 2:
        spo, spt, spth, spf, spfi = Lload_list


'''THE CHART'''

data = [[epo, ept, epth, epf, epfi],#8kg predicted length (e.g. epo = 'E'ight kg 'P'redicted test set 'O'ne)
        [eao, eat, eath, eaf, eafi],#8kg actual           (e.g. sath = 'S'ixteen kg 'A'ctual test set 'TH'ree)
        [spo, spt, spth, spf, spfi],#16kg predicted
        [sao, sat, sath, saf, safi]]#16kg actual

X = np.arange(5) # Five sets of bars (One per member of test set)
fig = plt.figure()
ax = fig.add_axes([0.1,0.2,0.85,0.7])
ax.bar(X + 0.0, data[0], color = 'b', width = 0.1, label = 'Precicted', zorder = 1) # 8kg length
ax.bar(X + 0.1, data[1], color = 'r', width = 0.1, label = 'Actual', zorder = 1) # 8kg length
ax.bar(X + 0.3, data[2], color = 'b', width = 0.1, zorder = 1) # 16kg length
ax.bar(X + 0.4, data[3], color = 'r', width = 0.1, zorder = 1) # 16kg length
ax.grid(which='major', axis='y')

plt.ylabel('Loaded Length (mm)', fontsize = Gconsts['FontSize'])
plt.title(Gconsts['Labels'][Gconsts['VertNum']], fontsize = Gconsts['FontSize']+4)
plt.ylim(6,20)

# Creates Green cross bars
x_pos = [0.05, 0.35, 1.05, 1.35, 2.05, 2.35, 3.05, 3.35, 4.05, 4.35]
ax.scatter(x_pos, Lo_Bar, color = 'green', marker = '_', s = 450, zorder = 2, label = 'Unloaded Length')

# x axis labels
bars = ['Test 1 8kg', 'Test 1 16kg', 'Test 2 8kg', 'Test 2 16kg', 'Test 3 8kg', 'Test 3 16kg',
        'Test 4 8kg', 'Test 4 16kg', 'Test 5 8kg', 'Test 5 16kg' ]
plt.xticks(x_pos, bars, rotation = 30,fontsize = Gconsts['FontSize']-2)

plt.legend(fontsize = Gconsts['FontSize']-2)

# Saving the plot
if Gconsts['Save'] == True:
    name = 'TestBars'+Gconsts['Labels'][Gconsts['VertNum']]+'.png'
    os.chdir('I:\MPhys Project\Research\Clean Dis References\Results\TestBars')
    plt.savefig(name, dpi = 200)
    os.chdir('I:\MPhys Project\Data\LoadedSpineData')